#include "cppdefs.h"
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
    real(rk)              :: dlon, dlat, dx, dy ! Could be arrays
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
    integer, intent(in)          :: nx, ny
    character(len=*), intent(in) :: topofile
    character(len=*), intent(in) :: lon, lat, bathy ! Variable names
    integer                      :: ii, jj

    dbghead(init_domain)

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
            d%seamask(ii, jj) = BEACH
          else
            d%seamask(ii, jj) = SEA
          end if
        else
          d%seamask(ii, jj) = LAND
        end if
      end do
    end do
    do ii = 1, nx
      if (d%depdata(ii, 1) .gt. ZERO) then
        d%seamask(ii, 1) = BOUNDARY
      else
        d%seamask(ii, 1) = LAND
      end if
      if (d%depdata(ii, ny) .gt. ZERO) then
        d%seamask(ii, ny) = BOUNDARY
      else
        d%seamask(ii, ny) = LAND
      end if
    end do
    do jj = 1, ny
      if (d%depdata(1, jj) .gt. ZERO) then
        d%seamask(1, jj) = BOUNDARY
      else
        d%seamask(1, jj) = LAND
      end if
      if (d%depdata(nx, jj) .gt. ZERO) then
        d%seamask(nx, jj) = BOUNDARY
      else
        d%seamask(nx, jj) = LAND
      end if
    end do

    FMT2, sum(d%seamask, mask=d%seamask == LAND), " land points"
    FMT2, sum(d%seamask, mask=d%seamask == SEA) / SEA, " sea points"
    FMT2, sum(d%seamask, mask=d%seamask == BEACH) / BEACH, " beach points"
    FMT2, sum(d%seamask, mask=d%seamask == BOUNDARY) / BOUNDARY, " boundary points"
    FMT2, nx * ny, " total points"

#ifdef DEBUG
#define FNAME "seamask.nc"
    DBG, "Saving seamask"
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

    dbgtail(init_domain)
    return
  end function ctor_domain
  !===========================================
  function get_lons_whole(this) result(res)
    class(t_domain), intent(in)  :: this
    real(rk), dimension(this%nx) :: res

    dbghead(get_lons_whole)

    res = this%lons

    dbgtail(get_lons_whole)
    return
  end function get_lons_whole
  !===========================================
  real(rk) function get_lons_idx(this, idx) result(res)
    class(t_domain), intent(in) :: this
    integer, intent(in) :: idx

    res = this%lons(idx)

    return
  end function get_lons_idx
  !===========================================
  real(rk) function get_lons_interp(this, idx) result(res)
    class(t_domain), intent(in) :: this
    real(rk), intent(in) :: idx
    real(rk) :: i0, i1

    dbghead(get_lons_interp)

    debug(idx)

    i0 = float(floor(idx)); debug(i0)
    i1 = float(floor(idx) + 1); debug(i1)

    debug(this%lons(int(i0)))
    debug(this%lons(int(i1)))

    res = (i1 - idx) / (i1 - i0) * this%lons(int(i0)) + (idx - i0) / (i1 - i0) * this%lons(int(i1))

    debug(res)

    dbgtail(get_lons_interp)
    return
  end function get_lons_interp
  !===========================================
  function get_lats_whole(this) result(res)
    class(t_domain), intent(in)  :: this
    real(rk), dimension(this%ny) :: res

    dbghead(get_lats_whole)

    res = this%lats

    dbgtail(get_lats_whole)
    return
  end function get_lats_whole
  !===========================================
  real(rk) function get_lats_idx(this, i) result(res)
    class(t_domain), intent(in) :: this
    integer, intent(in) :: i

    res = this%lats(i)

    return
  end function get_lats_idx
  !===========================================
  real(rk) function get_lats_interp(this, idx) result(res)
    class(t_domain), intent(in) :: this
    real(rk), intent(in) :: idx
    real(rk) :: i0, i1

    dbghead(get_lats_interp)

    debug(idx)

    i0 = float(floor(idx)); debug(i0)
    i1 = float(floor(idx) + 1); debug(i1)

    debug(this%lats(int(i0)))
    debug(this%lats(int(i1)))

    res = (i1 - idx) / (i1 - i0) * this%lats(int(i0)) + (idx - i0) / (i1 - i0) * this%lats(int(i1))

    debug(res)

    dbgtail(get_lats_interp)
    return
  end function get_lats_interp
  !===========================================
  function get_bathymetry_whole(this) result(res)
    class(t_domain), intent(in)           :: this
    real(rk), dimension(this%nx, this%ny) :: res

    res = this%depdata

    return
  end function get_bathymetry_whole
  !===========================================
  function get_bathymetry_idx(this, i, j) result(res)
    class(t_domain), intent(in) :: this
    integer, intent(in)         :: i, j
    real(rk)                    :: res

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
    class(t_domain), intent(in)          :: this
    integer, dimension(this%nx, this%ny) :: res

    res = this%seamask

  end function get_seamask_whole
  !===========================================
  function get_seamask_idx(this, i, j) result(res)
    class(t_domain), intent(in) :: this
    integer, intent(in)         :: i, j
    integer                     :: res

    res = this%seamask(i, j)

    return
  end function get_seamask_idx
  !===========================================
  subroutine lonlat2xy(this, lon, lat, x, y)
    class(t_domain), intent(in) :: this
    real(rk), intent(in)        :: lon, lat
    real(rk), intent(out)       :: x, y

    x = (lon - this%lon0) / this%dlon * this%dx
    y = (lat - this%lat0) / this%dlat * this%dy

    return
  end subroutine lonlat2xy
  !===========================================
  subroutine xy2lonlat(this, x, y, lon, lat)
    class(t_domain), intent(in) :: this
    real(rk), intent(in)        :: x, y
    real(rk), intent(out)       :: lon, lat

    lon = x / this%dx * this%dlon + this%lon0
    lat = y / this%dy * this%dlat + this%lat0

    return
  end subroutine xy2lonlat
  !===========================================
  subroutine get_indices_2d(this, lon, lat, i, j, ir, jr)
    !---------------------------------------------
    ! Should integer indices be int(irt) or nint(irt) (nearest)?
    !---------------------------------------------
    class(t_domain), intent(in)     :: this
    real(rk), intent(in)            :: lon, lat
    integer, optional, intent(out)  :: i, j
    real(rk), optional, intent(out) :: ir, jr
    real(rk)                        :: irt, jrt

    dbghead(get_indices_2d)

    debug(lon); debug(lat)

    irt = (lon - this%lon0) / this%dlon
    jrt = (lat - this%lat0) / this%dlat

#ifdef SNAP_TO_BOUNDS
    ! This is probably only necessary in case of large time steps
    if (irt <= ONE) irt = ONE
    if (irt > real(this%nx, rk)) irt = real(this%nx, rk)

    if (jrt <= ONE) jrt = ONE
    if (jrt > real(this%ny, rk)) jrt = real(this%ny, rk)
#endif

    debug(irt)
    debug(jrt)

    if (present(i)) then
      i = int(irt); debug(i)
      if (i < 1) call throw_error("domain :: get_indices_2d", "i is less than 1")
      if (i > this%nx) call throw_error("domain :: get_indices_2d", "i is greater than nx")
    end if

    if (present(j)) then
      j = int(jrt); debug(j)
      if (j < 1) call throw_error("domain :: get_indices_2d", "j is less than 1")
      if (j > this%ny) call throw_error("domain :: get_indices_2d", "j is greater than ny")
    end if

    if (present(ir)) then
      ir = irt; debug(ir)
      if (ir < ONE) call throw_error("domain :: get_indices_2d", "ir is less than 1")
      if (ir > real(this%nx, rk)) call throw_error("domain :: get_indices_2d", "ir is greater than nx")
    end if

    if (present(jr)) then
      jr = jrt; debug(jr)
      if (jr < ONE) call throw_error("domain :: get_indices_2d", "jr is less than 1")
      if (jr > real(this%ny, rk)) call throw_error("domain :: get_indices_2d", "jr is greater than ny")
    end if

    dbgtail(get_indices_2d)
    return
  end subroutine get_indices_2d

end module mod_domain
