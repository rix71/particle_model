#include "cppdefs.h"
module errors
  !----------------------------------------------------------------
  ! Basically stop
  !----------------------------------------------------------------
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: throw_error, throw_warning
  !===================================================
contains
  !===========================================
  subroutine throw_error(func, msg, code)

    character(len=*), intent(in)  :: func, msg
    integer, intent(in), optional :: code

    ERROR, trim(func)//": "//trim(msg)
    if (present(code)) then
      call perror("   ERROR: "//trim(func))
      ERROR, "err ", code
    end if
    stop

    return
  end subroutine throw_error
  !===========================================
  subroutine throw_warning(func, msg)

    character(len=*), intent(in) :: func, msg

    WARNING, trim(func)//": "//trim(msg)

    return
  end subroutine throw_warning

end module errors
