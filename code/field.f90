#include "cppdefs.h"
module mod_field
  !----------------------------------------------------------------
  ! This defines the field types that are used in the linked list.
  ! The field types are: static and dynamic. Both 1D, 2D, and 3D.
  ! The static field does not change in time (e.g. the bathymetry or depth axis)
  ! The dynamic field changes in time (e.g. the velocity or temperature field) and
  ! therefore defines two arrays: one for the current time step and one for the and one for the next time step.
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  use mod_variable, only: t_variable
  use mod_interp
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_field_static, t_field_dynamic
  public :: t_field_static_1d, t_field_static_2d, t_field_static_3d
  public :: t_field_dynamic_1d, t_field_dynamic_2d, t_field_dynamic_3d
  !---------------------------------------------
  type, abstract, extends(t_variable) :: t_field_static
  contains
    procedure :: print => print_info_static
    procedure, public :: gradient => gradient_static
    procedure, public :: slice => slice_static
  end type t_field_static
  !---------------------------------------------
  type, extends(t_field_static) :: t_field_static_1d
    private
    real(rk), dimension(:), allocatable :: data
    integer :: n
  contains
    private
    procedure, public :: init => init_1d_static
    procedure :: interpolate => interpolate_1d_static
    procedure, public :: set => set_value_1d_static
    procedure :: get_value_whole_1d_static, get_value_xy_1d_static, get_value_idx_1d_static
    generic, public :: get => get_value_whole_1d_static, get_value_xy_1d_static, get_value_idx_1d_static
  end type t_field_static_1d
  !---------------------------------------------
  type, extends(t_field_static) :: t_field_static_2d
    private
    real(rk), dimension(:, :), allocatable :: data
    integer :: n1, n2
  contains
    private
    procedure, public :: init => init_2d_static
    procedure :: interpolate => interpolate_2d_static
    procedure, public :: set => set_value_2d_static
    procedure :: get_value_whole_2d_static, get_value_xy_2d_static, get_value_idx_2d_static
    generic, public :: get => get_value_whole_2d_static, get_value_xy_2d_static, get_value_idx_2d_static
  end type t_field_static_2d
  !---------------------------------------------
  type, extends(t_field_static) :: t_field_static_3d
    private
    real(rk), dimension(:, :, :), allocatable :: data
    integer :: n1, n2, n3
    logical :: nan_bottom = .false., nan_top = .false.
  contains
    private
    procedure, public :: init => init_3d_static
    procedure :: interpolate => interpolate_3d_static
    procedure, public :: set => set_value_3d_static
    procedure :: get_value_whole_3d_static, get_value_xy_3d_static, get_value_idx_3d_static
    generic, public :: get => get_value_whole_3d_static, get_value_xy_3d_static, get_value_idx_3d_static
    procedure, public :: top_is_nan => top_is_nan_static
    procedure, public :: bottom_is_nan => bottom_is_nan_static
  end type t_field_static_3d
  !---------------------------------------------
  interface t_field_static
    module procedure ctor_field_static_1d
    module procedure ctor_field_static_2d
    module procedure ctor_field_static_3d
  end interface t_field_static
  !---------------------------------------------
  type, abstract, extends(t_variable) :: t_field_dynamic
    private
    logical :: set_t1 = .false.
    logical :: set_t2 = .false.
    real(rk) :: timestep
  contains
    private
    procedure, public :: init => init_dynamic
    procedure, public :: set_timestep, get_timestep
    procedure :: time_interpolation_field, time_interpolation_scalar, time_interpolation_1d, time_interpolation_2d, time_interpolation_3d
    generic :: time_interpolation => time_interpolation_field, time_interpolation_scalar, time_interpolation_1d, time_interpolation_2d, time_interpolation_3d
    procedure :: gradient => gradient_dynamic
    procedure, public :: slice => slice_dynamic
    procedure, public :: swap
  end type t_field_dynamic
  !---------------------------------------------
  type, extends(t_field_dynamic) :: t_field_dynamic_1d
    private
    type(t_field_static_1d), allocatable :: data_t1, data_t2
    integer :: n
  contains
    private
    procedure, public :: set => set_value_1d_dynamic
    procedure :: get_value_xy_1d_dynamic, get_value_idx_1d_dynamic
    generic, public :: get => get_value_xy_1d_dynamic, get_value_idx_1d_dynamic
  end type t_field_dynamic_1d
  !---------------------------------------------
  type, extends(t_field_dynamic) :: t_field_dynamic_2d
    private
    type(t_field_static_2d), allocatable :: data_t1, data_t2
    integer :: n1, n2
  contains
    private
    procedure, public :: set => set_value_2d_dynamic
    procedure :: get_value_xy_2d_dynamic, get_value_idx_2d_dynamic
    generic, public:: get => get_value_xy_2d_dynamic, get_value_idx_2d_dynamic
  end type t_field_dynamic_2d
  !---------------------------------------------
  type, extends(t_field_dynamic) :: t_field_dynamic_3d
    private
    type(t_field_static_3d), allocatable :: data_t1, data_t2
    integer :: n1, n2, n3
    logical :: nan_bottom = .false., nan_top = .false.
  contains
    private
    procedure, public :: set => set_value_3d_dynamic
    procedure :: get_value_xy_3d_dynamic, get_value_idx_3d_dynamic
    generic, public :: get => get_value_xy_3d_dynamic, get_value_idx_3d_dynamic
    procedure, public :: top_is_nan => top_is_nan_dynamic
    procedure, public :: bottom_is_nan => bottom_is_nan_dynamic
  end type t_field_dynamic_3d
  !---------------------------------------------
  interface t_field_dynamic
    module procedure ctor_field_dynamic_1d
    module procedure ctor_field_dynamic_2d
    module procedure ctor_field_dynamic_3d
  end interface t_field_dynamic

  !===================================================
contains
  !===========================================
  ! CONSTRUCTORS
  !===========================================
  type(t_field_static_1d) function ctor_field_static_1d(n, name)
    integer, intent(in) :: n
    character(len=*), intent(in) :: name

    call ctor_field_static_1d%init(n, name=name)

    return
  end function ctor_field_static_1d
  !===========================================
  type(t_field_static_2d) function ctor_field_static_2d(n1, n2, name)
    integer, intent(in) :: n1, n2
    character(len=*), intent(in) :: name

    call ctor_field_static_2d%init(n1, n2, name=name)

    return
  end function ctor_field_static_2d
  !===========================================
  type(t_field_static_3d) function ctor_field_static_3d(n1, n2, n3, name)
    integer, intent(in) :: n1, n2, n3
    character(len=*), intent(in) :: name

    call ctor_field_static_3d%init(n1, n2, n3, name=name)

    return
  end function ctor_field_static_3d
  !===========================================
  type(t_field_dynamic_1d) function ctor_field_dynamic_1d(n, timestep, name)
    integer, intent(in) :: n
    real(rk), intent(in) :: timestep
    character(len=*), intent(in) :: name

    call ctor_field_dynamic_1d%init(1, [n], name)
    call ctor_field_dynamic_1d%set_timestep(timestep)

    return
  end function ctor_field_dynamic_1d
  !===========================================
  type(t_field_dynamic_2d) function ctor_field_dynamic_2d(n1, n2, timestep, name)
    integer, intent(in) :: n1, n2
    real(rk), intent(in) :: timestep
    character(len=*), intent(in) :: name

    call ctor_field_dynamic_2d%init(2, [n1, n2], name)
    call ctor_field_dynamic_2d%set_timestep(timestep)

    return
  end function ctor_field_dynamic_2d
  !===========================================
  type(t_field_dynamic_3d) function ctor_field_dynamic_3d(n1, n2, n3, timestep, name)
    integer, intent(in) :: n1, n2, n3
    real(rk), intent(in) :: timestep
    character(len=*), intent(in) :: name

    dbghead(ctor_field_dynamic_3d)

    call ctor_field_dynamic_3d%init(3, [n1, n2, n3], name)
    call ctor_field_dynamic_3d%set_timestep(timestep)

    dbgtail(ctor_field_dynamic_3d)
    return
  end function ctor_field_dynamic_3d
  !===========================================
  ! INIT ROUTINES
  !===========================================
  subroutine init_dynamic(this, ndim, sdim, name, units)
    class(t_field_dynamic), intent(inout) :: this
    integer, intent(in) :: ndim
    integer, dimension(ndim), intent(in) :: sdim
    character(len=*), intent(in), optional :: name, units

    dbghead(field :: init_dynamic)

    select type (this)
    type is (t_field_dynamic_1d)
      DBG, "TYPE IS 1D DYNAMIC"
      if (ndim /= 1) call throw_error("field :: init_dynamic", "type of "//trim(name)//" is t_field_dynamic_1d but ndim /= 1")
      allocate (this%data_t1)
      allocate (this%data_t2)
      call this%data_t1%init(sdim(1))
      call this%data_t2%init(sdim(1))

      this%n = sdim(1)

      if (present(units)) then
        call this%data_t1%set_units(units)
        call this%data_t2%set_units(units)
      end if
      if (present(name)) then
        call this%data_t1%set_name(name)
        call this%data_t2%set_name(name)
      end if

    type is (t_field_dynamic_2d)
      DBG, "TYPE IS 2D DYNAMIC"
      if (ndim /= 2) call throw_error("field :: init_dynamic", "type of "//trim(name)//" is t_field_dynamic_2d but ndim /= 2")
      allocate (this%data_t1)
      allocate (this%data_t2)
      call this%data_t1%init(sdim(1), sdim(2))
      call this%data_t2%init(sdim(1), sdim(2))

      this%n1 = sdim(1)
      this%n2 = sdim(2)

      if (present(units)) then
        call this%data_t1%set_units(units)
        call this%data_t2%set_units(units)
      end if
      if (present(name)) then
        call this%data_t1%set_name(name)
        call this%data_t2%set_name(name)
      end if

    type is (t_field_dynamic_3d)
      DBG, "TYPE IS 3D DYNAMIC"
      if (ndim /= 3) call throw_error("field :: init_dynamic", "type of "//trim(name)//" is t_field_dynamic_3d but ndim /= 3")
      allocate (this%data_t1)
      allocate (this%data_t2)
      call this%data_t1%init(sdim(1), sdim(2), sdim(3))
      call this%data_t2%init(sdim(1), sdim(2), sdim(3))

      this%n1 = sdim(1)
      this%n2 = sdim(2)
      this%n3 = sdim(3)

      if (present(units)) then
        call this%data_t1%set_units(units)
        call this%data_t2%set_units(units)
      end if
      if (present(name)) then
        call this%data_t1%set_name(name)
        call this%data_t2%set_name(name)
      end if
    end select

    call this%set_dim(ndim)

    if (present(units)) then
      call this%set_units(units)
    end if

    if (present(name)) then
      call this%set_name(name)
    end if

    dbgtail(field :: init_dynamic)
    return
  end subroutine init_dynamic

  !===========================================
  subroutine init_1d_static(this, n, name, value, units)
    class(t_field_static_1d), intent(inout) :: this
    integer, intent(in) :: n
    real(rk), intent(in), optional :: value(n)
    character(len=*), intent(in), optional :: name, units

    allocate (this%data(n))
    if (present(value)) then
      this%data = value
    else
      this%data = ZERO
    end if
    this%n = n
    call this%set_dim(1)
    if (present(name)) then
      call this%set_name(name)
    end if
    if (present(units)) then
      call this%set_units(units)
    end if

    return
  end subroutine init_1d_static
  !===========================================
  subroutine init_2d_static(this, n1, n2, name, value, units)
    class(t_field_static_2d), intent(inout) :: this
    integer, intent(in) :: n1, n2
    real(rk), intent(in), optional :: value(n1, n2)
    character(len=*), intent(in), optional :: name, units

    allocate (this%data(n1, n2))
    if (present(value)) then
      this%data = value
    else
      this%data = ZERO
    end if
    this%n1 = n1
    this%n2 = n2
    call this%set_dim(2)
    if (present(name)) then
      call this%set_name(name)
    end if
    if (present(units)) then
      call this%set_units(units)
    end if

    return
  end subroutine init_2d_static
  !===========================================
  subroutine init_3d_static(this, n1, n2, n3, name, value, units)
    class(t_field_static_3d), intent(inout) :: this
    integer, intent(in) :: n1, n2, n3
    real(rk), intent(in), optional :: value(n1, n2, n3)
    character(len=*), intent(in), optional :: name, units

    dbghead(field :: init_3d_static)

    DBG, "ALLOCATING"
    allocate (this%data(n1, n2, n3))
    DBG, "ALLOCATED"
    if (present(value)) then
      DBG, "PRESENT VALUE"
      this%data = value
    else
      DBG, "NOT PRESENT VALUE"
      this%data = ZERO
    end if
    this%n1 = n1
    this%n2 = n2
    this%n3 = n3
    DBG, "SET DIM"
    call this%set_dim(3)
    if (present(name)) then
      DBG, "PRESENT NAME"
      call this%set_name(name)
    end if
    if (present(units)) then
      DBG, "PRESENT UNITS"
      call this%set_units(units)
    end if

    dbgtail(field :: init_3d_static)
    return
  end subroutine init_3d_static
  !===========================================
  ! Interpolation routines
  !===========================================
  real(rk) function interpolate_1d_static(this, x) result(res)
    class(t_field_static_1d), intent(in) :: this
    real(rk), intent(in) :: x ! The index to interpolate at
    real(rk) :: x1, x2
    integer :: i1, i2

    i1 = floor(x)
    i2 = i1 + 1
    x1 = real(i1, rk)
    x2 = real(i2, rk)

    call linearinterp(x1, x2, this%data(i1), this%data(i2), x, res)

    return
  end function interpolate_1d_static
  !===========================================
  real(rk) function interpolate_2d_static(this, x, y) result(res)
    class(t_field_static_2d), intent(in) :: this
    real(rk), intent(in) :: x, y ! The indices to interpolate at
    real(rk) :: c11, c12, c21, c22
    integer :: i, j
    real(rk) :: x1, x2, y1, y2

    i = floor(x)
    j = floor(y)

    x1 = real(i, rk)
    x2 = real(i + 1, rk)
    y1 = real(j, rk)
    y2 = real(j + 1, rk)

    c11 = this%data(i, j)
    c12 = this%data(i, j + 1)
    c21 = this%data(i + 1, j)
    c22 = this%data(i + 1, j + 1)

    call bilinearinterp(x1, x1, x2, x2, y1, y2, c11, c12, c21, c22, x, y, res)

    return
  end function interpolate_2d_static
  !===========================================
  real(rk) function interpolate_3d_static(this, x, y, z) result(res)
    class(t_field_static_3d), intent(in) :: this
    real(rk), intent(in) :: x, y, z ! The indices to interpolate at
    real(rk) :: c111, c121, c211, c221, c112, c122, c212, c222
    integer :: i1, j1, k1, i2, j2, k2
    real(rk) :: x1, x2, y1, y2, z1, z2
    logical :: do_bilin

    dbghead(field :: interpolate_3d_static)

    do_bilin = .false.

    i1 = floor(x)
    j1 = floor(y)
    k1 = floor(z)

    if (i1 >= this%n1 .or. j1 >= this%n2 .or. &
        i1 <= 1 .or. j1 <= 1) then
      res = ZERO
      return
    end if

    if (k1 == this%n3) then
      if (this%n3 == 1) then
        ! It's a 3D field, but only one layer in the vertical, so we'll just do bilinear interpolation
        do_bilin = .true.
      else
        k1 = k1 - 1
      end if
    else if (k1 < 1) then
      k1 = 1
    end if

    i2 = i1 + 1
    j2 = j1 + 1
    k2 = k1 + 1

    x1 = real(i1, rk)
    x2 = real(i2, rk)
    y1 = real(j1, rk)
    y2 = real(j2, rk)
    z1 = real(k1, rk)
    z2 = real(k2, rk)

    c111 = this%data(i1, j1, k1); debug(c111)
    c121 = this%data(i1, j2, k1); debug(c121)
    c211 = this%data(i2, j1, k1); debug(c211)
    c221 = this%data(i2, j2, k1); debug(c221)

    if (do_bilin) then
      call bilinearinterp(x1, x1, x2, x2, y1, y2, c111, c121, c211, c221, x, y, res)
      debug(res)
      dbgtail(field :: interpolate_3d_static)
      return
    end if

    c112 = this%data(i1, j1, k2)
    c122 = this%data(i1, j2, k2)
    c212 = this%data(i2, j1, k2)
    c222 = this%data(i2, j2, k2)

    call trilinearinterp(x1, x2, y1, y2, z1, z2, c111, c121, c211, c221, c112, c122, c212, c222, x, y, z, res)
    debug(res)

    dbgtail(field :: interpolate_3d_static)
    return
  end function interpolate_3d_static
  !===========================================
  function time_interpolation_field(this, t, f1, f2) result(res)
    class(t_field_dynamic), intent(in) :: this
    real(rk), intent(in) :: t
    class(t_field_static), intent(in) :: f1, f2
    class(t_field_static), allocatable :: res
    real(rk), dimension(:), allocatable :: c11, c12
    real(rk), dimension(:, :), allocatable :: c21, c22
    real(rk), dimension(:, :, :), allocatable :: c31, c32

    if (f1%get_dim() /= f2%get_dim()) then
      call throw_error("field :: time_interpolation_field", "The fields to interpolate between must have the same dimension")
    end if

    select type (f1)
    type is (t_field_static_1d)
      call f1%get(c11)
    type is (t_field_static_2d)
      call f1%get(c21)
    type is (t_field_static_3d)
      call f1%get(c31)
    end select

    select type (f2)
    type is (t_field_static_1d)
      call f2%get(c12)
    type is (t_field_static_2d)
      call f2%get(c22)
    type is (t_field_static_3d)
      call f2%get(c32)
    end select

    select type (f1)
    type is (t_field_static_1d)
      allocate (res, source=f1)
      select type (res)
      type is (t_field_static_1d)
        call res%set(c11 + (c12 - c11) * (t / this%timestep))
      end select
    type is (t_field_static_2d)
      allocate (res, source=f1)
      select type (res)
      type is (t_field_static_2d)
        call res%set(c21 + (c22 - c21) * (t / this%timestep))
      end select
    type is (t_field_static_3d)
      allocate (res, source=f1)
      select type (res)
      type is (t_field_static_3d)
        call res%set(c31 + (c32 - c31) * (t / this%timestep))
      end select
    end select

    return
  end function time_interpolation_field
  !===========================================
  real(rk) function time_interpolation_scalar(this, t, f1, f2) result(res)
    class(t_field_dynamic), intent(in) :: this
    real(rk), intent(in) :: t
    real(rk), intent(in) :: f1, f2

    res = f1 + (f2 - f1) * (t / this%timestep)

    return
  end function time_interpolation_scalar
  !===========================================
  function time_interpolation_1d(this, t, f1, f2, n) result(res)
    class(t_field_dynamic), intent(in) :: this
    real(rk), intent(in) :: t
    real(rk), dimension(:), intent(in) :: f1, f2
    integer, intent(in) :: n
    real(rk), dimension(n) :: res

    res = f1 + (f2 - f1) * (t / this%timestep)

    return
  end function time_interpolation_1d
  !===========================================
  function time_interpolation_2d(this, t, f1, f2, n1, n2) result(res)
    class(t_field_dynamic), intent(in) :: this
    real(rk), intent(in) :: t
    real(rk), dimension(:, :), intent(in) :: f1, f2
    integer, intent(in) :: n1, n2
    real(rk), dimension(n1, n2) :: res

    res = f1 + (f2 - f1) * (t / this%timestep)

    return
  end function time_interpolation_2d
  !===========================================
  function time_interpolation_3d(this, t, f1, f2, n1, n2, n3) result(res)
    class(t_field_dynamic), intent(in) :: this
    real(rk), intent(in) :: t
    real(rk), dimension(:, :, :), intent(in) :: f1, f2
    integer, intent(in) :: n1, n2, n3
    real(rk), dimension(n1, n2, n3) :: res

    res = f1 + (f2 - f1) * (t / this%timestep)

    return
  end function time_interpolation_3d
  !===========================================
  ! GETTER FUNCTIONS
  !===========================================
  subroutine get_value_whole_1d_static(this, res)
    class(t_field_static_1d), intent(in) :: this
    real(rk), dimension(:), intent(out) :: res

    res = this%data

    return
  end subroutine get_value_whole_1d_static
  !===========================================
  subroutine get_value_whole_2d_static(this, res)
    class(t_field_static_2d), intent(in) :: this
    real(rk), dimension(:, :), intent(out) :: res

    res = this%data

    return
  end subroutine get_value_whole_2d_static
  !===========================================
  subroutine get_value_whole_3d_static(this, res)
    class(t_field_static_3d), intent(in) :: this
    real(rk), dimension(:, :, :), intent(out) :: res

    res = this%data

    return
  end subroutine get_value_whole_3d_static
  !===========================================
  subroutine get_value_xy_1d_static(this, x, res)
    class(t_field_static_1d), intent(in) :: this
    real(rk), intent(in) :: x
    real(rk), intent(out) :: res

    res = this%interpolate(x)

    return
  end subroutine get_value_xy_1d_static
  !===========================================
  subroutine get_value_xy_2d_static(this, x, y, res)
    class(t_field_static_2d), intent(in) :: this
    real(rk), intent(in) :: x, y
    real(rk), intent(out) :: res

    res = this%interpolate(x, y)

    return
  end subroutine get_value_xy_2d_static
  !===========================================
  subroutine get_value_xy_3d_static(this, x, y, z, res)
    class(t_field_static_3d), intent(in) :: this
    real(rk), intent(in) :: x, y, z
    real(rk), intent(out) :: res

    res = this%interpolate(x, y, z)

    return
  end subroutine get_value_xy_3d_static
  !===========================================
  subroutine get_value_idx_1d_static(this, i, res)
    class(t_field_static_1d), intent(in) :: this
    integer, intent(in) :: i
    real(rk), intent(out) :: res

    res = this%data(i)

    return
  end subroutine get_value_idx_1d_static
  !===========================================
  subroutine get_value_idx_2d_static(this, i, j, res)
    class(t_field_static_2d), intent(in) :: this
    integer, intent(in) :: i, j
    real(rk), intent(out) :: res

    res = this%data(i, j)

    return
  end subroutine get_value_idx_2d_static
  !===========================================
  subroutine get_value_idx_3d_static(this, i, j, k, res)
    class(t_field_static_3d), intent(in) :: this
    integer, intent(in) :: i, j, k
    real(rk), intent(out) :: res

    res = this%data(i, j, k)

    return
  end subroutine get_value_idx_3d_static
  !===========================================
  subroutine get_value_xy_1d_dynamic(this, t, x, res)
    class(t_field_dynamic_1d), intent(in) :: this
    real(rk), intent(in) :: t, x
    real(rk), intent(out) :: res
    real(rk) :: c1, c2

    call this%data_t1%get(x, c1)
    call this%data_t2%get(x, c2)

    res = this%time_interpolation(t, c1, c2)

    return
  end subroutine get_value_xy_1d_dynamic
  !===========================================
  subroutine get_value_xy_2d_dynamic(this, t, x, y, res)
    class(t_field_dynamic_2d), intent(in) :: this
    real(rk), intent(in) :: t, x, y
    real(rk), intent(out) :: res
    real(rk) :: c1, c2

    call this%data_t1%get(x, y, c1)
    call this%data_t2%get(x, y, c2)

    res = this%time_interpolation(t, c1, c2)

    return
  end subroutine get_value_xy_2d_dynamic
  !===========================================
  subroutine get_value_xy_3d_dynamic(this, t, x, y, z, res)
    class(t_field_dynamic_3d), intent(in) :: this
    real(rk), intent(in) :: t, x, y, z
    real(rk), intent(out) :: res
    real(rk) :: c1, c2

    call this%data_t1%get(x, y, z, c1)
    call this%data_t2%get(x, y, z, c2)

    res = this%time_interpolation(t, c1, c2)

    return
  end subroutine get_value_xy_3d_dynamic
  !===========================================
  subroutine get_value_idx_1d_dynamic(this, t, i, res)
    class(t_field_dynamic_1d), intent(in) :: this
    real(rk), intent(in) :: t
    integer, intent(in) :: i
    real(rk), intent(out) :: res
    real(rk) :: c1, c2

    call this%data_t1%get(i, c1)
    call this%data_t2%get(i, c2)

    res = this%time_interpolation(t, c1, c2)

    return
  end subroutine get_value_idx_1d_dynamic
  !===========================================
  subroutine get_value_idx_2d_dynamic(this, t, i, j, res)
    class(t_field_dynamic_2d), intent(in) :: this
    real(rk), intent(in) :: t
    integer, intent(in) :: i, j
    real(rk), intent(out) :: res
    real(rk) :: c1, c2

    call this%data_t1%get(i, j, c1)
    call this%data_t2%get(i, j, c2)

    res = this%time_interpolation(t, c1, c2)

    return
  end subroutine get_value_idx_2d_dynamic
  !===========================================
  subroutine get_value_idx_3d_dynamic(this, t, i, j, k, res)
    class(t_field_dynamic_3d), intent(in) :: this
    real(rk), intent(in) :: t
    integer, intent(in) :: i, j, k
    real(rk), intent(out) :: res
    real(rk) :: c1, c2

    call this%data_t1%get(i, j, k, c1)
    call this%data_t2%get(i, j, k, c2)

    res = this%time_interpolation(t, c1, c2)

    return
  end subroutine get_value_idx_3d_dynamic
  !===========================================
  real(rk) function get_timestep(this) result(res)
    class(t_field_dynamic), intent(in) :: this

    res = this%timestep

    return
  end function get_timestep
  !===========================================
  ! SETTER FUNCTIONS
  !===========================================
  subroutine set_value_1d_static(this, field)
    class(t_field_static_1d), intent(inout) :: this
    real(rk), intent(in) :: field(:)

    this%data = field

    return
  end subroutine set_value_1d_static
  !===========================================
  subroutine set_value_2d_static(this, field)
    class(t_field_static_2d), intent(inout) :: this
    real(rk), intent(in) :: field(:, :)

    this%data = field

    return
  end subroutine set_value_2d_static
  !===========================================
  subroutine set_value_3d_static(this, field, t_nan, b_nan)
    class(t_field_static_3d), intent(inout) :: this
    real(rk), intent(in) :: field(:, :, :)
    logical, intent(in), optional :: t_nan, b_nan

    this%data = field

    if (present(t_nan)) this%nan_top = t_nan
    if (present(b_nan)) this%nan_bottom = b_nan

    return
  end subroutine set_value_3d_static
  !===========================================
  subroutine set_value_1d_dynamic(this, field)
    class(t_field_dynamic_1d), intent(inout) :: this
    real(rk), intent(in) :: field(:)

    if (.not. this%set_t1) then
      call this%data_t1%set(field)
      this%set_t1 = .true.
    else if (.not. this%set_t2) then
      call this%data_t2%set(field)
      this%set_t2 = .true.
    else
      call throw_error("field :: set_value_1d_dynamic", "Both t1 and t2 are already set.")
    end if

    return
  end subroutine set_value_1d_dynamic
  !===========================================
  subroutine set_value_2d_dynamic(this, field)
    class(t_field_dynamic_2d), intent(inout) :: this
    real(rk), intent(in) :: field(:, :)

    if (.not. this%set_t1) then
      call this%data_t1%set(field)
      this%set_t1 = .true.
    else if (.not. this%set_t2) then
      call this%data_t2%set(field)
      this%set_t2 = .true.
    else
      call throw_error("field :: set_value_2d_dynamic", "Both t1 and t2 are already set.")
    end if

    return
  end subroutine set_value_2d_dynamic
  !===========================================
  subroutine set_value_3d_dynamic(this, field, t_nan, b_nan)
    class(t_field_dynamic_3d), intent(inout) :: this
    real(rk), intent(in) :: field(:, :, :)
    logical, intent(in), optional :: t_nan, b_nan

    if (.not. this%set_t1) then
      call this%data_t1%set(field)
      this%set_t1 = .true.
    else if (.not. this%set_t2) then
      call this%data_t2%set(field)
      this%set_t2 = .true.
    else
      call throw_error("field :: set_value_3d_dynamic", "Both t1 and t2 are already set.")
    end if

    if (present(t_nan)) this%nan_top = t_nan
    if (present(b_nan)) this%nan_bottom = b_nan

    return
  end subroutine set_value_3d_dynamic
  !===========================================
  subroutine set_timestep(this, timestep)
    class(t_field_dynamic), intent(inout) :: this
    real(rk), intent(in) :: timestep

    this%timestep = timestep

    return
  end subroutine set_timestep
  !===========================================
  ! SWAP FUNCTIONS
  !===========================================
  subroutine swap(this)
    class(t_field_dynamic), intent(inout) :: this

    if (this%set_t1 .and. this%set_t2) then
      select type (this)
      type is (t_field_dynamic_1d)
        this%data_t1 = this%data_t2
      type is (t_field_dynamic_2d)
        this%data_t1 = this%data_t2
      type is (t_field_dynamic_3d)
        this%data_t1 = this%data_t2
      end select
      this%set_t2 = .false.
    end if

    return
  end subroutine swap
  !===========================================
  ! SLICE FUNCTIONS
  !===========================================
  subroutine slice_static(this, slice_dim, idx_other, res)
    !---------------------------------------------
    ! Returns a slice of the field in the 'slice_dim' dimension (1D field).
    ! The 'idx_other' array contains the indices of the other dimensions.
    !---------------------------------------------
    class(t_field_static), intent(in) :: this
    integer, intent(in) :: slice_dim
    integer, intent(in) :: idx_other(:)
    type(t_field_static_1d), intent(out) :: res

    select type (this)
    type is (t_field_static_1d)
      call throw_error("field :: slice_static", "Cannot slice a 1D field.")
    type is (t_field_static_2d)
      if (size(idx_other) /= 1) then
        call throw_error("field :: slice_static", "Invalid number of indices.")
      end if
      select case (slice_dim)
      case (1)
        call res%init(n=this%n1, &
                      name="slice_"//trim(this%get_name()), &
                      value=this%data(:, idx_other(1)))
      case (2)
        call res%init(n=this%n2, &
                      name="slice_"//trim(this%get_name()), &
                      value=this%data(idx_other(1), :))
      case default
        call throw_error("field :: slice_static", "Invalid dimension.")
      end select
    type is (t_field_static_3d)
      if (size(idx_other) /= 2) then
        call throw_error("field :: slice_static", "Invalid number of indices.")
      end if
      select case (slice_dim)
      case (1)
        call res%init(n=this%n1, &
                      name="slice_"//trim(this%get_name()), &
                      value=this%data(:, idx_other(1), idx_other(2)))
      case (2)
        call res%init(n=this%n2, &
                      name="slice_"//trim(this%get_name()), &
                      value=this%data(idx_other(1), :, idx_other(2)))
      case (3)
        call res%init(n=this%n3, &
                      name="slice_"//trim(this%get_name()), &
                      value=this%data(idx_other(1), idx_other(2), :))
      case default
        call throw_error("field :: slice_static", "Invalid dimension.")
      end select
    end select
    return
  end subroutine slice_static
  !===========================================
  subroutine slice_dynamic(this, slice_dim, t, idx_other, res)
    class(t_field_dynamic), intent(in) :: this
    integer, intent(in) :: slice_dim
    real(rk), intent(in) :: t
    integer, intent(in) :: idx_other(:)
    type(t_field_static_1d), intent(out) :: res
    type(t_field_static_1d) :: f1, f2
    class(t_field_static), allocatable :: f_interp

    select type (this)
    type is (t_field_dynamic_1d)
      call throw_error("field :: slice_dynamic", "Cannot slice a 1D field.")

    type is (t_field_dynamic_2d)
      call this%data_t1%slice(slice_dim, idx_other, f1)
      call this%data_t2%slice(slice_dim, idx_other, f2)
      f_interp = this%time_interpolation(t, f1, f2)
      select type (f_interp)
      type is (t_field_static_1d)
        res = f_interp
      end select

    type is (t_field_dynamic_3d)
      call this%data_t1%slice(slice_dim, idx_other, f1)
      call this%data_t2%slice(slice_dim, idx_other, f2)
      f_interp = this%time_interpolation(t, f1, f2)
      select type (f_interp)
      type is (t_field_static_1d)
        res = f_interp
      end select
    end select

    return
  end subroutine slice_dynamic
  !===========================================
  ! GRADIENT FUNCTIONS
  !===========================================
  subroutine gradient_static(this, dim, res)
    class(t_field_static), intent(in) :: this
    integer, intent(in) :: dim
    class(t_field_static), allocatable, intent(out) :: res

    select type (this)
    type is (t_field_static_1d)
      allocate (t_field_static_1d :: res)
      select type (res)
      type is (t_field_static_1d)
        call res%init(n=this%n - 1, &
                      name="gradient_"//trim(this%get_name()), &
                      value=this%data(2:this%n) - this%data(1:this%n - 1))
      end select

    type is (t_field_static_2d)
      allocate (t_field_static_2d :: res)
      select type (res)
      type is (t_field_static_2d)
        select case (dim)
        case (1)
          call res%init(n1=this%n1 - 1, n2=this%n2, &
                        name="gradient_"//trim(this%get_name()), &
                        value=this%data(2:this%n1, :) - this%data(1:this%n1 - 1, :))
        case (2)
          call res%init(n1=this%n1, n2=this%n2 - 1, &
                        name="gradient_"//trim(this%get_name()), &
                        value=this%data(:, 2:this%n2) - this%data(:, 1:this%n2 - 1))
        case default
          call throw_error("field :: gradient", "dim must be 1 or 2.")
        end select
      end select

    type is (t_field_static_3d)
      allocate (t_field_static_3d :: res)
      select type (res)
      type is (t_field_static_3d)
        select case (dim)
        case (1)
          call res%init(n1=this%n1 - 1, n2=this%n2, n3=this%n3, &
                        name="gradient_"//trim(this%get_name()), &
                        value=this%data(2:this%n1, :, :) - this%data(1:this%n1 - 1, :, :))
        case (2)
          call res%init(n1=this%n1, n2=this%n2 - 1, n3=this%n3, &
                        name="gradient_"//trim(this%get_name()), &
                        value=this%data(:, 2:this%n2, :) - this%data(:, 1:this%n2 - 1, :))
        case (3)
          call res%init(n1=this%n1, n2=this%n2, n3=this%n3 - 1, &
                        name="gradient_"//trim(this%get_name()), &
                        value=this%data(:, :, 2:this%n3) - this%data(:, :, 1:this%n3 - 1))
        case default
          call throw_error("field :: gradient", "dim must be 1, 2 or 3.")
        end select
      end select
    end select

    return
  end subroutine gradient_static
  !===========================================
  subroutine gradient_dynamic(this, dim, res)
    class(t_field_dynamic), intent(in) :: this
    integer, intent(in) :: dim
    class(t_field_dynamic), allocatable, intent(out) :: res

    call throw_error("field :: gradient_dynamic", "Not implemented yet.")

    return
  end subroutine gradient_dynamic
  !===========================================
  ! INFO FUNCTIONS (no idea how else to call these)
  !===========================================
  logical function top_is_nan_static(this)
    class(t_field_static_3d), intent(in) :: this

    top_is_nan_static = this%nan_top

    return
  end function top_is_nan_static
  !===========================================
  logical function bottom_is_nan_static(this)
    class(t_field_static_3d), intent(in) :: this

    bottom_is_nan_static = this%nan_bottom

    return
  end function bottom_is_nan_static
  !===========================================
  logical function top_is_nan_dynamic(this)
    class(t_field_dynamic_3d), intent(in) :: this

    top_is_nan_dynamic = this%nan_top

    return
  end function top_is_nan_dynamic
  !===========================================
  logical function bottom_is_nan_dynamic(this)
    class(t_field_dynamic_3d), intent(in) :: this

    bottom_is_nan_dynamic = this%nan_bottom

    return
  end function bottom_is_nan_dynamic
  !===========================================
  subroutine print_info_static(this)
    class(t_field_static), intent(in) :: this

    FMT3, "Static field:"
    call this%print_metadata()
    select type (this)
    type is (t_field_static_1d)
      FMT3, "n = ", this%n
    type is (t_field_static_2d)
      FMT3, "n1 = ", this%n1, " n2 = ", this%n2
    type is (t_field_static_3d)
      FMT3, "n1 = ", this%n1, " n2 = ", this%n2, " n3 = ", this%n3
    end select

  end subroutine print_info_static
  !===========================================

end module mod_field
