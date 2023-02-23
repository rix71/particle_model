#include "cppdefs.h"
module mod_variable
  !----------------------------------------------------------------
  ! This is a base class for all field data types.
  !----------------------------------------------------------------
  use mod_precdefs
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_variable
  !---------------------------------------------
  type, abstract :: t_variable
    private
    character(len=LEN_CHAR_S) :: name = ""
    character(len=LEN_CHAR_S) :: unit = ""
    integer :: dim = 0
    real(rk) :: missing_value = -9999.0_rk
  contains
    procedure :: set_metadata
    procedure :: get_name, get_units, get_dim
    procedure :: set_name, set_units, set_dim
    procedure :: set_missing_value, get_missing_value
    procedure :: print_metadata
  end type t_variable
  !===================================================
contains
  !===========================================
  subroutine set_name(this, name)
    class(t_variable), intent(inout) :: this
    character(len=*), intent(in) :: name

    this%name = name

    return
  end subroutine set_name
  !===========================================
  subroutine set_units(this, unit)
    class(t_variable), intent(inout) :: this
    character(len=*), intent(in) :: unit

    this%unit = unit

    return
  end subroutine set_units
  !===========================================
  subroutine set_dim(this, dim)
    class(t_variable), intent(inout) :: this
    integer, intent(in) :: dim

    this%dim = dim

    return
  end subroutine set_dim
  !===========================================
  subroutine set_missing_value(this, missing_value)
    class(t_variable), intent(inout) :: this
    real(rk), intent(in) :: missing_value

    this%missing_value = missing_value

    return
  end subroutine set_missing_value
  !===========================================
  function get_name(this) result(name)
    class(t_variable), intent(in) :: this
    character(len=LEN_CHAR_S) :: name

    name = this%name

    return
  end function get_name
  !===========================================
  function get_units(this) result(unit)
    class(t_variable), intent(in) :: this
    character(len=LEN_CHAR_S) :: unit

    unit = this%unit

    return
  end function get_units
  !===========================================
  function get_dim(this) result(dim)
    class(t_variable), intent(in) :: this
    integer :: dim

    dim = this%dim

    return
  end function get_dim
  !===========================================
  function get_missing_value(this) result(missing_value)
    class(t_variable), intent(in) :: this
    real(rk) :: missing_value

    missing_value = this%missing_value

    return
  end function get_missing_value
  !===========================================
  subroutine set_metadata(this, name, unit, dim, missing_value)
    class(t_variable), intent(inout) :: this
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: unit
    integer, intent(in), optional :: dim
    real(rk), intent(in), optional :: missing_value

    if (present(name)) then
      call this%set_name(name)
    end if
    if (present(unit)) then
      call this%set_units(unit)
    end if
    if (present(dim)) then
      call this%set_dim(dim)
    end if
    if (present(missing_value)) then
      call this%set_missing_value(missing_value)
    end if

    return
  end subroutine set_metadata
  !===========================================
  subroutine print_metadata(this)
    class(t_variable), intent(in) :: this

    FMT2, "Variable metadata:"
    FMT3, "Name: "//trim(this%get_name())
    FMT3, "Units: "//trim(this%get_units())
    FMT3, "Dim: ", this%get_dim()
    FMT3, "Missing value: ", this%get_missing_value()

    return
  end subroutine print_metadata
end module mod_variable
