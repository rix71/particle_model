#if (defined(SNAP_TO_BOUNDS) && defined(DEBUG))
#warning SNAP_TO_BOUNDS defined, indices can be modified
#endif
#include "cppdefs.h"
module mod_field
  use mod_errors
  use mod_precdefs
  use mod_interp, only: bilinearinterp, trilinearinterp
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_field
  !---------------------------------------------
  type t_field
    private
    real(rk), allocatable     :: data_t1(:, :, :)
    real(rk), allocatable     :: data_t2(:, :, :)
    integer                   :: nx, ny, nz
    integer, public           :: n_dims
    character(len=LEN_CHAR_S) :: nc_varname = "field"
    logical                   :: set_t1 = .false.
    logical                   :: set_t2 = .false.
    real(rk)                  :: timestep

  contains
    private
    generic, public   :: get => get_interp, get_nointerp
    procedure         :: get_interp, get_nointerp
    ! procedure, public :: get => get_interp, get_nointerp
    procedure, public :: set
    procedure, public :: swap
    procedure, public :: get_varname
    procedure, public :: get_array_pointer
    procedure, public :: gradient
    generic           :: time_interp => time_interp_scalar, time_interp_vector, time_interp_matrix
    procedure         :: time_interp_scalar, time_interp_vector, time_interp_matrix
    procedure, public :: get_array_1D
    final             :: dtor_field
  end type t_field
  !---------------------------------------------
  interface t_field
    module procedure :: ctor_field
  end interface t_field
  !---------------------------------------------
  integer :: ierr
  !===================================================
contains
  !===========================================
  type(t_field) function ctor_field(nx, ny, nz, timestep, nc_varname) result(g)

    integer, intent(in)                    :: nx, ny
    integer, optional, intent(in)          :: nz
    real(rk), intent(in)                   :: timestep
    character(len=*), optional, intent(in) :: nc_varname

    if (present(nz)) then
      g%nx = nx
      g%ny = ny
      g%nz = nz
      g%n_dims = 3
      allocate (g%data_t1(nx, ny, nz), stat=ierr)
      if (ierr .ne. 0) call throw_error("field :: field", "Could not allocate", ierr)
      allocate (g%data_t2(nx, ny, nz), stat=ierr)
      if (ierr .ne. 0) call throw_error("field :: field", "Could not allocate", ierr)
    else
      g%nx = nx
      g%ny = ny
      g%nz = 1
      g%n_dims = 2
      allocate (g%data_t1(nx, ny, 1), stat=ierr)
      if (ierr .ne. 0) call throw_error("field :: field", "Could not allocate", ierr)
      allocate (g%data_t2(nx, ny, 1), stat=ierr)
      if (ierr .ne. 0) call throw_error("field :: field", "Could not allocate", ierr)
    end if
    g%data_t1 = ZERO
    g%data_t2 = ZERO
    g%timestep = timestep
    if (present(nc_varname)) g%nc_varname = nc_varname

  end function ctor_field
  !===========================================
  subroutine dtor_field(this)
    type(t_field), intent(inout) :: this

    if (allocated(this%data_t1)) deallocate (this%data_t1, stat=ierr)
    if (ierr .ne. 0) call throw_warning("field :: field destructor", "Could not deallocate data (t1)")
    if (allocated(this%data_t2)) deallocate (this%data_t2, stat=ierr)
    if (ierr .ne. 0) call throw_warning("field :: field destructor", "Could not deallocate data (t2)")

  end subroutine dtor_field
  !===========================================
  character(len=LEN_CHAR_S) function get_varname(this) result(res)
    class(t_field) :: this

    res = this%nc_varname

  end function get_varname
  !===========================================
  real(rk) function time_interp_scalar(this, f1, f2, t) result(res)
    class(t_field), intent(in) :: this
    real(rk), intent(in) :: f1, f2, t

    res = f1 + (f2 - f1) * (t / this%timestep)

    return
  end function time_interp_scalar
  !===========================================
  function time_interp_vector(this, f1, f2, t, n) result(res)
    class(t_field), intent(in) :: this
    integer, intent(in) :: n
    real(rk), dimension(n), intent(in) :: f1, f2
    real(rk), intent(in) :: t
    real(rk) :: res(n)

    res = f1 + (f2 - f1) * (t / this%timestep)

    return
  end function time_interp_vector
  !===========================================
  function time_interp_matrix(this, f1, f2, t, n, m, k) result(res)
    class(t_field), intent(in) :: this
    integer, intent(in) :: n, m, k
    real(rk), dimension(n, m, k), intent(in) :: f1, f2
    real(rk), intent(in) :: t
    real(rk) :: res(n, m, k)

    res = f1 + (f2 - f1) * (t / this%timestep)

    return
  end function time_interp_matrix
  !===========================================
  subroutine get_interp(this, t, x, y, z, &
                        res)
    !---------------------------------------------
    ! Getter function
    ! Get (interpolated) values from the field using real indices
    ! Input t: time between t1 and t2 (t1 <= t <= t2 or 0 <= t <= timestep)
    !---------------------------------------------
    class(t_field), intent(in)     :: this
    real(rk), intent(in)           :: t
    real(rk) AIM_INTENT            :: x, y
    real(rk), optional, intent(in) :: z
    ! integer                        :: inbr, iadd, jadd
    real(rk), intent(out)          :: res
    integer                        :: i, j, k
    real(rk)                       :: f1, f2
    real(rk)                       :: x1, x2, &
                                      y1, y2, &
                                      z1, z2, &
                                      c11, c12, c21, c22, &
                                      c111, c121, c211, c221, c112, c122, c212, c222

    dbghead(get_interp)

    debug(t); debug(x); debug(y)

    i = floor(x); debug(i)
    x1 = float(floor(x)); debug(x1)
    x2 = float(floor(x) + 1); debug(x2)
    if (x >= real(this%nx, rk)) then
#ifdef SNAP_TO_BOUNDS
      ! Right boundary
      ! can only be trusted when the SNAP_TO_BOUNDS flag is defined
      DBG, "Snapping x down"
      i = this%nx - 1; debug(i)
      x = real(this%nx, rk); debug(x)
      x1 = real(this%nx, rk) - ONE; debug(x1)
      x2 = real(this%nx, rk); debug(x2)
#else
      call throw_error("field :: get_interp", "Index (i) out of bounds!")
#endif
    end if
    if (x < ONE) then
#ifdef SNAP_TO_BOUNDS
      ! Left boundary
      ! can only be trusted when the SNAP_TO_BOUNDS flag is defined
      DBG, "Snapping x up"
      i = 1; debug(i)
      x = ONE; debug(x)
      x1 = ONE; debug(x1)
      x2 = 2.*ONE; debug(x2)
#else
      call throw_error("field :: get_interp", "Index (i) out of bounds!")
#endif
    end if

    j = floor(y); debug(j)
    y1 = float(floor(y)); debug(y1)
    y2 = float(floor(y) + 1); debug(y2)
    if (y >= real(this%ny, rk)) then
#ifdef SNAP_TO_BOUNDS
      ! Top boundary
      DBG, "Snapping y down"
      j = this%ny - 1; debug(j)
      y = real(this%ny, rk); debug(y)
      y1 = real(this%ny, rk) - ONE; debug(y1)
      y2 = real(this%ny, rk); debug(y2)
#else
      call throw_error("field :: get_interp", "Index (j) out of bounds!")
#endif
    end if
    if (y < ONE) then
#ifdef SNAP_TO_BOUNDS
      ! Bottom boundary
      DBG, "Snapping y up"
      j = 1; debug(j)
      y = ONE; debug(y)
      y1 = ONE; debug(y1)
      y2 = 2.*ONE; debug(y2)
#else
      call throw_error("field :: get_interp", "Index (j) out of bounds!")
#endif
    end if

    if (present(z)) then
      debug(z)
      DBG, "Spatial interpolation 3D"

      k = floor(z); debug(k)
      z1 = float(floor(z)); debug(z1)
      z2 = float(floor(z) + 1); debug(z2)

      if (k == this%nz) then
        ! Use the lower layer if the particle is exactly at the surface
        k = k - 1
        z1 = z1 - ONE
        z2 = z2 - ONE
      end if

      c111 = this%data_t1(i, j, k); debug(c111)
      c121 = this%data_t1(i, j + 1, k); debug(c121)
      c211 = this%data_t1(i + 1, j, k); debug(c211)
      c221 = this%data_t1(i + 1, j + 1, k); debug(c221)
      c112 = this%data_t1(i, j, k + 1); debug(c112)
      c122 = this%data_t1(i, j + 1, k + 1); debug(c122)
      c212 = this%data_t1(i + 1, j, k + 1); debug(c212)
      c222 = this%data_t1(i + 1, j + 1, k + 1); debug(c222)

      call trilinearinterp(x1, x2, y1, y2, z1, z2, c111, c121, c211, c221, c112, c122, c212, c222, x, y, z, f1)
      debug(f1)

      c111 = this%data_t2(i, j, k); debug(c111)
      c121 = this%data_t2(i, j + 1, k); debug(c121)
      c211 = this%data_t2(i + 1, j, k); debug(c211)
      c221 = this%data_t2(i + 1, j + 1, k); debug(c221)
      c112 = this%data_t2(i, j, k + 1); debug(c112)
      c122 = this%data_t2(i, j + 1, k + 1); debug(c122)
      c212 = this%data_t2(i + 1, j, k + 1); debug(c212)
      c222 = this%data_t2(i + 1, j + 1, k + 1); debug(c222)

      call trilinearinterp(x1, x2, y1, y2, z1, z2, c111, c121, c211, c221, c112, c122, c212, c222, x, y, z, f2)
      debug(f2)
    else
      DBG, "Spatial interpolation 2D"

      c11 = this%data_t1(i, j, 1); debug(c11)
      c12 = this%data_t1(i, j + 1, 1); debug(c12)
      c21 = this%data_t1(i + 1, j, 1); debug(c21)
      c22 = this%data_t1(i + 1, j + 1, 1); debug(c22)

      call bilinearinterp(x1, x1, x2, x2, y1, y2, c11, c12, c21, c22, x, y, f1)
      debug(f1)

      c11 = this%data_t2(i, j, 1)
      c12 = this%data_t2(i, j + 1, 1)
      c21 = this%data_t2(i + 1, j, 1)
      c22 = this%data_t2(i + 1, j + 1, 1)

      call bilinearinterp(x1, x1, x2, x2, y1, y2, c11, c12, c21, c22, x, y, f2)
      debug(f2)
    end if

    res = this%time_interp(f1, f2, t)
    debug(res)

    dbgtail(get_interp)
    return
  end subroutine get_interp
  !===========================================
  subroutine get_nointerp(this, t, i, j, k, res)
    !---------------------------------------------
    ! Getter function with integer indices
    !---------------------------------------------
    class(t_field), intent(in) :: this
    real(rk), intent(in) :: t
    integer, intent(in) :: i, j
    integer, optional, intent(in) :: k
    real(rk), intent(out) :: res
    real(rk) :: f1, f2

    dbghead(get_nointerp)

    if (present(k)) then
      f1 = this%data_t1(i, j, k)
      f2 = this%data_t2(i, j, k)
    else
      f1 = this%data_t1(i, j, 1)
      f2 = this%data_t2(i, j, 1)
    end if
    res = this%time_interp(f1, f2, t)

    dbgtail(get_nointerp)
    return
  end subroutine get_nointerp
  !===========================================
  subroutine get_array_1D(this, t, i, j, n, dim, res)
    class(t_field), intent(in) :: this
    real(rk), intent(in) :: t
    integer, intent(in) :: i, j
    integer, intent(in) :: n
    integer, optional, intent(in) :: dim
    real(rk), intent(out) :: res(n)
    real(rk) :: f1(n), f2(n)
    integer :: idim

    dbghead(get_array_1D)

    debug(t); debug(i); debug(j); debug(n)

    if (present(dim)) then
      debug(dim)
      idim = dim
    else
      idim = 1
    end if

    select case (idim)
    case (1)
      if (n .ne. size(this%data_t1, dim=1)) call throw_error("field :: get_array_1D", "n does not match the size of dimension 1")
      f1 = this%data_t1(:, i, j)
      f2 = this%data_t2(:, i, j)
    case (2)
      if (n .ne. size(this%data_t1, dim=2)) call throw_error("field :: get_array_1D", "n does not match the size of dimension 2")
      f1 = this%data_t1(i, :, j)
      f2 = this%data_t2(i, :, j)
    case (3)
      if (n .ne. size(this%data_t1, dim=3)) call throw_error("field :: get_array_1D", "n does not match the size of dimension 3")
      f1 = this%data_t1(i, j, :)
      f2 = this%data_t2(i, j, :)
    case default
      call throw_error("field :: get_array_1D", "Dimension must be 1, 2 or 3!")
    end select

    res = this%time_interp(f1, f2, t, n)

    dbgtail(get_array_1D)
    return
  end subroutine get_array_1D
  !===========================================
  subroutine set(this, data)
    class(t_field), intent(inout) :: this
    real(rk), intent(in) :: data(this%nx, this%ny, this%nz)

    if (.not. this%set_t1) then
      DBG, "--> T1"
      this%data_t1 = data
      this%set_t1 = .true.
    else if (.not. this%set_t2) then
      DBG, "--> T2"
      this%data_t2 = data
      this%set_t2 = .true.
    else
      call throw_error("field :: set", "Fields T1 and T2 already exist!")
    end if

    return
  end subroutine set
  !===========================================
  function get_array_pointer(this) result(p_data)
    class(t_field), intent(inout) :: this
    real(rk), pointer :: p_data(:, :, :)

    if (.not. this%set_t1) then
      DBG, "--> T1"
      p_data = this%data_t1
      this%set_t1 = .true.
    else if (.not. this%set_t2) then
      DBG, "--> T2"
      p_data = this%data_t2
      this%set_t2 = .true.
    else
      call throw_error("field :: get_array_pointer", "Fields T1 and T2 already exist!")
    end if

    return
  end function get_array_pointer
  !===========================================
  function gradient(this, t, dx, dim) result(grad)
    class(t_field), intent(in) :: this
    real(rk), intent(in) :: t, dx
    integer, intent(in) :: dim
    real(rk), dimension(this%nx, this%ny, this%nz) :: grad
    real(rk) :: data_t(this%nx, this%ny, this%nz)

    data_t = this%time_interp(this%data_t1, this%data_t2, t, this%nx, this%ny, this%nz)

    select case (dim)
    case (1)
      grad(:this%nx - 1, :, :) = (data_t(2:, :, :) - data_t(:this%nx - 1, :, :))
      grad(this%nx, :, :) = grad(this%nx - 1, :, :)
    case (2)
      grad(:, :this%ny - 1, :) = (data_t(:, 2:, :) - data_t(:, :this%ny - 1, :))
      grad(:, this%ny, :) = grad(:, this%ny - 1, :)
    case default
      call throw_error("field :: gradient", "Dimension must be 1 or 2")
    end select

    grad = grad / dx

  end function gradient
  !===========================================
  subroutine swap(this)
    class(t_field), intent(inout) :: this

    dbghead(swap)

    if (this%set_t1 .and. this%set_t2) then
      DBG, trim(this%nc_varname)//": T2 --> T1"
      this%data_t1 = this%data_t2
      this%set_t2 = .false.
    end if

    dbgtail(swap)
    return
  end subroutine swap

end module mod_field
