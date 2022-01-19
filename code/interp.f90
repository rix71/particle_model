#include "cppdefs.h"
module interp
  !----------------------------------------------------------------
  ! This module contains different interpolation methods
  !----------------------------------------------------------------
  use precdefs
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: bilinearinterp, trilinearinterp, timeinterp
  !===================================================
contains
  !===========================================
  subroutine bilinearinterp(x11, x12, x21, x22, y1, y2, z11, z12, z21, z22, x, y, z)

    real, intent(in)      :: x11, x12, x21, x22
    real, intent(in)      :: y1, y2
    real(rk), intent(in)  :: z11, z12, z21, z22
    real(rk), intent(in)  :: x, y
    real(rk), intent(out) :: z
    real(rk)              :: z1, z2

    z1 = (x12 - x) / (x12 - x11) * z11 + (x - x11) / (x12 - x11) * z12
    z2 = (x22 - x) / (x22 - x21) * z21 + (x - x21) / (x22 - x21) * z22
    z = (y2 - y) / (y2 - y1) * z1 + (y - y1) / (y2 - y1) * z2

    return
  end subroutine bilinearinterp
  !===========================================
  subroutine trilinearinterp(x1, x2, y1, y2, z1, z2, &
                             c111, c121, c211, c221, c112, c122, c212, c222, &
                             x, y, z, c)
    !---------------------------------------------
    ! indexing: c_ijk
    ! e.g. 111 - front left bottom corner
    !      121 - front right bottom corner
    !      221 - back right bottom corner
    !      222 - back right top corner
    !       11 - front bottom corner along x-axis
    !       22 - back top corner along x-axis
    !        1 - bottom point along y-axis
    !        2 - top point along y-axis
    ! TODO: Maybe calculate dx, dy and dz for each edge?
    !---------------------------------------------

    real, intent(in)      :: x1, x2
    real, intent(in)      :: y1, y2
    real, intent(in)      :: z1, z2
    real(rk)              :: dx, dy, dz
    real(rk), intent(in)  :: c111, c121, c211, c221, c112, c122, c212, c222
    real(rk), intent(in)  :: x, y, z
    real(rk), intent(out) :: c
    real(rk)              :: c11, c12, c21, c22
    real(rk)              :: c1, c2

    dx = (x - x1) / (x2 - x1)
    dy = (y - y1) / (y2 - y1)
    dz = (z - z1) / (z2 - z1)

    c11 = c111 * (1 - dx) + c121 * dx
    c12 = c112 * (1 - dx) + c122 * dx
    c21 = c211 * (1 - dx) + c221 * dx
    c22 = c212 * (1 - dx) + c222 * dx

    c1 = c11 * (1 - dy) + c21 * dy
    c2 = c12 * (1 - dy) + c22 * dy

    c = c1 * (1 - dz) + c2 * dz

    return
  end subroutine trilinearinterp
  !===========================================
  subroutine timeinterp(arr_t1, arr_t2, arr_out, dt, dt_arr, nx, ny, nz, method)
    integer, intent(in)   :: nx, ny, nz
    integer, intent(in)   :: method
    real(rk), intent(in)  :: arr_t1(nx, ny, nz)
    real(rk), intent(in)  :: arr_t2(nx, ny, nz)
    real(rk), intent(in)  :: dt, dt_arr
    real(rk), intent(out) :: arr_out(nx, ny, nz)

    select case (method)
    case (0)
      arr_out = arr_t1
    case (1)
      arr_out = arr_t1 + (arr_t2 - arr_t1) / dt_arr * dt
    end select

  end subroutine timeinterp
end module interp
