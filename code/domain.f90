#include "cppdefs.h"
#include "field.h"
module mod_domain
  !----------------------------------------------------------------
  ! Initialise domain and seamask
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  use mod_interp, only: bilinearinterp
!   use mod_params, only: pi
!   use mod_domain_vars
  use nc_manager, only: nc_read_real_1d, nc_read_real_2d
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_domain
  !---------------------------------------------
  real(rk), parameter :: pi = 4.*atan(1.)
  !---------------------------------------------
  type t_domain
    private
    integer, public       :: nx, ny
    real(rk), allocatable :: lons(:), lats(:)
    real(rk), allocatable :: depdata(:, :)
    integer, allocatable  :: seamask(:, :)
    real(rk)              :: lon0, lat0, lon1, lat1
    real(rk), public      :: dlon, dlat, dx, dy ! Could be arrays
  contains
    private
    procedure, public :: lonlat2xy, xy2lonlat
    procedure, public :: get_indices_2d
    generic, public   :: get_bathymetry => get_bathymetry_whole, get_bathymetry_idx, get_bathymetry_idx_interp
    procedure         :: get_bathymetry_whole, get_bathymetry_idx, get_bathymetry_idx_interp
    generic, public   :: get_seamask => get_seamask_whole, get_seamask_idx
    procedure         :: get_seamask_whole, get_seamask_idx
    generic, public   :: get_lons => get_lons_whole, get_lons_idx, get_lons_interp
    procedure         :: get_lons_whole, get_lons_idx, get_lons_interp
    generic, public   :: get_lats => get_lats_whole, get_lats_idx, get_lats_interp
    procedure         :: get_lats_whole, get_lats_idx, get_lats_interp

  end type t_domain
  !---------------------------------------------
  interface t_domain
    module procedure :: ctor_domain
  end interface t_domain
  !===================================================
contains
  !===========================================
  type(t_domain) function ctor_domain(nx, ny, topofile, lon, lat, bathy) result(d)
    !---------------------------------------------
    ! Initialize the global longitude/latitude
    ! and seamask
    !---------------------------------------------
#ifdef DEBUG
    use nc_manager
    integer :: nc_x_dimid, nc_y_dimid
#endif
    integer, intent(in) :: nx, ny
    character(len=*), intent(in) :: topofile
    character(len=*), intent(in) :: lon, lat, bathy ! Variable names
    integer :: ii, jj

    d%nx = nx
    d%ny = ny

    allocate (d%lons(nx), d%lats(ny))

    FMT2, "Reading coordinates"
    call nc_read_real_1d(trim(TOPOFILE), trim(lon), nx, d%lons)
    call nc_read_real_1d(trim(TOPOFILE), trim(lat), ny, d%lats)

    d%lat0 = d%lats(1); d%lat1 = d%lats(ny); d%dlat = d%lats(2) - d%lats(1); d%dy = d%dlat * 60.*1852.
    d%lon0 = d%lons(1); d%lon1 = d%lons(nx); d%dlon = d%lons(2) - d%lons(1); d%dx = d%dlon * 60.*1852.*cos(0.5 * (d%lat0 + d%lat1) * pi / 180.)

    FMT2, LINE
    FMT2, "Coordinates:"
    FMT3, var2val(d%lat0), "[deg N], ", var2val(d%lat1), "[deg N]"
    FMT3, var2val(d%lon0), "[deg E], ", var2val(d%lon1), "[deg E]"
    FMT2, "Cell size:"
    FMT3, var2val(d%dlat), "[deg], ", var2val(d%dy), "[m]"
    FMT3, var2val(d%dlon), "[deg], ", var2val(d%dx), "[m]"

    allocate (d%depdata(nx, ny), d%seamask(nx, ny))

    FMT2, "Reading bathymetry"
    call nc_read_real_2d(trim(TOPOFILE), trim(bathy), nx, ny, d%depdata)

    !---------------------------------------------
    ! TODO: Seamask could have another value (4) to represent boundaries.
    !       Boundary should have a thickness!
    FMT2, "Making seamask"
    do ii = 2, nx - 1
      do jj = 2, ny - 1
        if (d%depdata(ii, jj) .gt. ZERO) then
          if ((d%depdata(ii + 1, jj) .le. 0.0) .or. (d%depdata(ii - 1, jj) .le. 0.0) .or. &
              (d%depdata(ii, jj + 1) .le. 0.0) .or. (d%depdata(ii, jj - 1) .le. 0.0) .or. &
              (d%depdata(ii + 1, jj + 1) .le. 0.0) .or. (d%depdata(ii + 1, jj - 1) .le. 0.0) .or. &
              (d%depdata(ii - 1, jj - 1) .le. 0.0) .or. (d%depdata(ii - 1, jj + 1) .le. 0.0)) then
            d%seamask(ii, jj) = DOM_BEACH
          else
            d%seamask(ii, jj) = DOM_SEA
          end if
        else
          d%seamask(ii, jj) = DOM_LAND
        end if
      end do
    end do
    do ii = 1, nx
      if (d%depdata(ii, 1) .gt. ZERO) then
        d%seamask(ii, 1) = DOM_BOUNDARY
      else
        d%seamask(ii, 1) = DOM_LAND
      end if
      if (d%depdata(ii, ny) .gt. ZERO) then
        d%seamask(ii, ny) = DOM_BOUNDARY
      else
        d%seamask(ii, ny) = DOM_LAND
      end if
    end do
    do jj = 1, ny
      if (d%depdata(1, jj) .gt. ZERO) then
        d%seamask(1, jj) = DOM_BOUNDARY
      else
        d%seamask(1, jj) = DOM_LAND
      end if
      if (d%depdata(nx, jj) .gt. ZERO) then
        d%seamask(nx, jj) = DOM_BOUNDARY
      else
        d%seamask(nx, jj) = DOM_LAND
      end if
    end do

    FMT2, sum(d%seamask, mask=d%seamask == DOM_LAND), " land points"
    FMT2, sum(d%seamask, mask=d%seamask == DOM_SEA) / DOM_SEA, " sea points"
    FMT2, sum(d%seamask, mask=d%seamask == DOM_BEACH) / DOM_BEACH, " beach points"
    FMT2, sum(d%seamask, mask=d%seamask == DOM_BOUNDARY) / DOM_BOUNDARY, " boundary points"
    FMT2, nx * ny, " total points"

#ifdef DEBUG
#define FNAME "seamask.nc"
    call nc_initialise(FNAME)
    call nc_add_dimension(FNAME, "lon", nc_x_dimid, nx)
    call nc_add_dimension(FNAME, "lat", nc_y_dimid, ny)
    call nc_add_variable(FNAME, "seamask", "int", 2, [nc_x_dimid, nc_y_dimid])
    call nc_add_variable(FNAME, "lon", "float", 1, [nc_x_dimid])
    call nc_add_variable(FNAME, "lat", "float", 1, [nc_y_dimid])
    call nc_write(FNAME, d%seamask, "seamask", nx, ny)
    call nc_write(FNAME, d%lons, "lon", nx)
    call nc_write(FNAME, d%lats, "lat", ny)
#undef FNAME
#endif

    FMT2, "Finished init_domain"

    return
  end function ctor_domain
  !===========================================
  function get_lons_whole(this) result(res)
    class(t_domain), intent(in) :: this
    real(rk), dimension(this%nx):: res

    res = this%lons

    return
  end function get_lons_whole
  !===========================================
  real(rk) function get_lons_idx(this, idx) result(res)
    class(t_domain), intent(in) :: this
    integer AIM_INTENT :: idx

    if (idx < 1) then
#ifdef SNAP_TO_BOUNDS
      idx = 1
#else
      call throw_error("domain :: get_lons_idx", "Index out of bounds! (Less than 1)")
#endif
    end if
    if (idx > this%nx) then
#ifdef SNAP_TO_BOUNDS
      idx = this%nx
#else
      call throw_error("domain :: get_lons_idx", "Index out of bounds! (Greater than nx)")
#endif
    end if

    res = this%lons(idx)

    return
  end function get_lons_idx
  !===========================================
  real(rk) function get_lons_interp(this, idx) result(res)
    class(t_domain), intent(in) :: this
    real(rk) AIM_INTENT :: idx
    real(rk) :: i0, i1

#ifdef SNAP_TO_BOUNDS
    if (idx < ONE) idx = ONE
    if (idx > real(this%nx, rk)) idx = real(this%nx, rk)
#else
    if (idx >= real(this%nx, rk)) then
      call throw_error("domain :: get_lons_interp", "Index out of bounds!")
    end if
#endif

    i0 = float(floor(idx)); 
    i1 = float(floor(idx) + 1); 
    if (int(i1) <= this%nx) then
      res = (i1 - idx) / (i1 - i0) * this%lons(int(i0)) + (idx - i0) / (i1 - i0) * this%lons(int(i1))
    else
      ! I guess this is a kinda acceptable solution to index errors
      res = (i1 - idx) / (i1 - i0) * this%lons(int(i0)) + (idx - i0) / (i1 - i0) * (this%lons(int(i0)) + this%dlon)
    end if

    return
  end function get_lons_interp
  !===========================================
  function get_lats_whole(this) result(res)
    class(t_domain), intent(in)  :: this
    real(rk), dimension(this%ny) :: res

    res = this%lats

    return
  end function get_lats_whole
  !===========================================
  real(rk) function get_lats_idx(this, idx) result(res)
    class(t_domain), intent(in) :: this
    integer AIM_INTENT :: idx

    if (idx < 1) then
#ifdef SNAP_TO_BOUNDS
      idx = 1
#else
      call throw_error("domain :: get_lats_idx", "Index out of bounds! (Less than 1)")
#endif
    end if
    if (idx > this%nx) then
#ifdef SNAP_TO_BOUNDS
      idx = this%nx
#else
      call throw_error("domain :: get_lats_idx", "Index out of bounds! (Greater than ny)")
#endif
    end if

    res = this%lats(idx)

    return
  end function get_lats_idx
  !===========================================
  real(rk) function get_lats_interp(this, idx) result(res)
    class(t_domain), intent(in) :: this
    real(rk) AIM_INTENT :: idx
    real(rk) :: i0, i1

#ifdef SNAP_TO_BOUNDS
    if (idx < ONE) idx = ONE
    if (idx > real(this%ny, rk)) idx = real(this%ny, rk)
#else
    if ((idx >= real(this%ny, rk)) .or. (idx < ONE)) then
      call throw_error("domain :: get_lats_interp", "Index out of bounds!")
    end if
#endif

    i0 = float(floor(idx)); 
    i1 = float(floor(idx) + 1); 
    if (int(i1) <= this%ny) then
      res = (i1 - idx) / (i1 - i0) * this%lats(int(i0)) + (idx - i0) / (i1 - i0) * this%lats(int(i1))
    else
      res = (i1 - idx) / (i1 - i0) * this%lats(int(i0)) + (idx - i0) / (i1 - i0) * (this%lats(int(i0)) + this%dlat)
    end if

    return
  end function get_lats_interp
  !===========================================
  function get_bathymetry_whole(this) result(res)
    class(t_domain), intent(in) :: this
    real(rk), dimension(this%nx, this%ny) :: res

    res = this%depdata

    return
  end function get_bathymetry_whole
  !===========================================
  function get_bathymetry_idx(this, i, j) result(res)
    class(t_domain), intent(in) :: this
    integer AIM_INTENT :: i, j
    real(rk) :: res

#ifdef SNAP_TO_BOUNDS
    if (i < 1) i = 1
    if (i > this%nx) i = this%nx
    if (j < 1) j = 1
    if (j > this%ny) j = this%ny
#endif

    if (i < 1) call throw_error("domain :: get_bathymetry_idx", "i is less than 1")
    if (i > this%nx) call throw_error("domain :: get_bathymetry_idx", "i is greater than nx")
    if (j < 1) call throw_error("domain :: get_bathymetry_idx", "j is less than 1")
    if (j > this%ny) call throw_error("domain :: get_bathymetry_idx", "j is greater than ny")

    res = this%depdata(i, j)

    return
  end function get_bathymetry_idx
  !===========================================
  function get_bathymetry_idx_interp(this, x, y) result(res)
    class(t_domain), intent(in) :: this
    real(rk), intent(in)        :: x, y ! Indices, not coordinates!
    real(rk)                    :: x1, x2, &
                                   y1, y2, &
                                   c11, c12, c21, c22
    integer                     :: i, j
    real(rk)                    :: res

    i = floor(x)
    x1 = float(floor(x))
    x2 = float(floor(x) + 1)

    j = floor(y)
    y1 = float(floor(y))
    y2 = float(floor(y) + 1)

    c11 = this%depdata(i, j)
    c12 = this%depdata(i, j + 1)
    c21 = this%depdata(i + 1, j)
    c22 = this%depdata(i + 1, j + 1)

    call bilinearinterp(x1, x1, x2, x2, y1, y2, c11, c12, c21, c22, x, y, res)

    return
  end function get_bathymetry_idx_interp
  !===========================================
  function get_seamask_whole(this) result(res)
    class(t_domain), intent(in) :: this
    integer, dimension(this%nx, this%ny) :: res

    res = this%seamask

  end function get_seamask_whole
  !===========================================
  function get_seamask_idx(this, i, j) result(res)
    class(t_domain), intent(in) :: this
    integer AIM_INTENT          :: i, j
    integer                     :: res

#ifdef SNAP_TO_BOUNDS
    if (i < 1) i = 1
    if (i > this%nx) i = this%nx
    if (j < 1) j = 1
    if (j > this%ny) j = this%ny
#endif

    if (i < 1) call throw_error("domain :: get_seamask_idx", "i is less than 1")
    if (i > this%nx) call throw_error("domain :: get_seamask_idx", "i is greater than nx")
    if (j < 1) call throw_error("domain :: get_seamask_idx", "j is less than 1")
    if (j > this%ny) call throw_error("domain :: get_seamask_idx", "j is greater than ny")

    res = this%seamask(i, j)

    return
  end function get_seamask_idx
  !===========================================
  subroutine lonlat2xy(this, lon, lat, x, y)
    class(t_domain), intent(in) :: this
    real(rk), intent(in)  :: lon, lat
    real(rk), intent(out) :: x, y

    x = (lon - this%lon0) / this%dlon * this%dx
    y = (lat - this%lat0) / this%dlat * this%dy

    return
  end subroutine lonlat2xy
  !===========================================
  subroutine xy2lonlat(this, x, y, lon, lat)
    class(t_domain), intent(in) :: this
    real(rk), intent(in) :: x, y
    real(rk), intent(out) :: lon, lat

    lon = x / this%dx * this%dlon + this%lon0
    lat = y / this%dy * this%dlat + this%lat0

    return
  end subroutine xy2lonlat
  !===========================================
  subroutine get_indices_2d(this, lon, lat, i, j, ir, jr)
    !---------------------------------------------
    ! Should integer indices be int(irt) or nint(irt) (nearest)?
    !---------------------------------------------
    class(t_domain), intent(in) :: this

    real(rk), intent(in)            :: lon, lat
    integer, optional, intent(out)  :: i, j
    real(rk), optional, intent(out) :: ir, jr
    integer                         :: it, jt
    real(rk)                        :: irt, jrt

    it = minloc(abs(this%lons - lon), dim=1); 
    jt = minloc(abs(this%lats - lat), dim=1); 
    ! if (this%lons(it) > lon) it = it - 1
    ! if (this%lats(jt) > lat) jt = jt - 1

    if (it == this%nx) then

      if (lon > this%lons(this%nx)) then
#ifdef SNAP_TO_BOUNDS

        irt = real(this%nx, rk)
#else
        call throw_error("domain :: get_indices_2d", "it is greater than nx")
#endif
      else
        if (this%lons(it) > lon) it = it - 1 ! This is probably always true
        irt = it + (lon - this%lons(it)) / (this%lons(it + 1) - this%lons(it))
      end if
    else if (it == 1) then

      if (lon < this%lons(1)) then
#ifdef SNAP_TO_BOUNDS

        irt = ONE
#else
        call throw_error("domain :: get_indices_2d", "it is less than 1")
#endif
      else
        ! if (this%lons(it) > lon) it = it - 1
        irt = it + (lon - this%lons(it)) / (this%lons(it + 1) - this%lons(it))
      end if
    else
      if (this%lons(it) > lon) it = it - 1
      irt = it + (lon - this%lons(it)) / (this%lons(it + 1) - this%lons(it))
    end if

    if (jt == this%ny) then

      if (lat > this%lats(this%ny)) then
#ifdef SNAP_TO_BOUNDS

        jrt = real(this%ny, rk)
#else
        call throw_error("domain :: get_indices_2d", "jt is greater than ny")
#endif
      else
        if (this%lats(jt) > lat) jt = jt - 1 ! Also always true
        jrt = jt + (lat - this%lats(jt)) / (this%lats(jt + 1) - this%lats(jt))
      end if
    else if (jt == 1) then

      if (lat < this%lats(1)) then
#ifdef SNAP_TO_BOUNDS

        jrt = ONE
#else
        call throw_error("domain :: get_indices_2d", "jt is less than 1")
#endif
      else
        ! if (this%lats(jt) > lat) jt = jt - 1
        jrt = jt + (lat - this%lats(jt)) / (this%lats(jt + 1) - this%lats(jt))
      end if
    else
      if (this%lats(jt) > lat) jt = jt - 1
      jrt = jt + (lat - this%lats(jt)) / (this%lats(jt + 1) - this%lats(jt))
    end if

    if (present(i)) then
      i = it; 
      if (i < 1) call throw_error("domain :: get_indices_2d", "i is less than 1")
      if (i > this%nx) call throw_error("domain :: get_indices_2d", "i is greater than nx")
    end if

    if (present(j)) then
      j = jt; 
      if (j < 1) call throw_error("domain :: get_indices_2d", "j is less than 1")
      if (j > this%ny) call throw_error("domain :: get_indices_2d", "j is greater than ny")
    end if

    if (present(ir)) then
      ir = irt; 
      if (ir < ONE) call throw_error("domain :: get_indices_2d", "ir is less than 1")
      if (ir > real(this%nx, rk)) call throw_error("domain :: get_indices_2d", "ir is greater than nx")
    end if

    if (present(jr)) then
      jr = jrt; 
      if (jr < ONE) call throw_error("domain :: get_indices_2d", "jr is less than 1")
      if (jr > real(this%ny, rk)) call throw_error("domain :: get_indices_2d", "jr is greater than ny")
    end if

    return
  end subroutine get_indices_2d
  !===========================================
!   subroutine get_indices_2d(this, lon, lat, i, j, ir, jr)
!     !---------------------------------------------
!     ! Should integer indices be int(irt) or nint(irt) (nearest)?
!     !---------------------------------------------
!     class(t_domain), intent(in) :: this

!     real(rk), intent(in)            :: lon, lat
!     integer, optional, intent(out)  :: i, j
!     real(rk), optional, intent(out) :: ir, jr
!     real(rk)                        :: irt, jrt

!

!

!     irt = (lon - this%lon0) / this%dlon + 1
!     jrt = (lat - this%lat0) / this%dlat + 1

! #ifdef SNAP_TO_BOUNDS
!     ! This is probably only necessary in case of large time steps
!     if (irt <= ONE) irt = ONE
!     if (irt > real(this%nx, rk)) irt = real(this%nx, rk)

!     if (jrt <= ONE) jrt = ONE
!     if (jrt > real(this%ny, rk)) jrt = real(this%ny, rk)
! #endif

!
!

!     if (present(i)) then
!       i = int(irt);
!       if (i < 1) call throw_error("domain :: get_indices_2d", "i is less than 1")
!       if (i > this%nx) call throw_error("domain :: get_indices_2d", "i is greater than nx")
!     end if

!     if (present(j)) then
!       j = int(jrt);
!       if (j < 1) call throw_error("domain :: get_indices_2d", "j is less than 1")
!       if (j > this%ny) call throw_error("domain :: get_indices_2d", "j is greater than ny")
!     end if

!     if (present(ir)) then
!       ir = irt;
!       if (ir < ONE) call throw_error("domain :: get_indices_2d", "ir is less than 1")
!       if (ir > real(this%nx, rk)) call throw_error("domain :: get_indices_2d", "ir is greater than nx")
!     end if

!     if (present(jr)) then
!       jr = jrt;
!       if (jr < ONE) call throw_error("domain :: get_indices_2d", "jr is less than 1")
!       if (jr > real(this%ny, rk)) call throw_error("domain :: get_indices_2d", "jr is greater than ny")
!     end if

!
!     return
!   end subroutine get_indices_2d

end module mod_domain
