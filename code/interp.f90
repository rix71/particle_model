#include "cppdefs.h"
module mod_interp
  !----------------------------------------------------------------
  ! This module contains different interpolation methods
  !----------------------------------------------------------------
  use mod_precdefs
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: bilinearinterp, trilinearinterp, linearinterp, nbrs
  !---------------------------------------------
  integer, parameter :: nbrs(2, 8) = reshape([0, 1, 1, 0, -1, 0, 0, -1, 1, 1, 1, -1, -1, -1, -1, 1], [2, 8])
  !===================================================
contains
  !===========================================
  subroutine linearinterp(x1, x2, y1, y2, x, y)
    real(rk), intent(in)  :: x1, x2
    real(rk), intent(in)  :: y1, y2
    real(rk), intent(in)  :: x
    real(rk), intent(out) :: y

    y = (y2 - y1) / (x2 - x1) * (x - x1) + y1

    return
  end subroutine linearinterp
  !===========================================
  subroutine bilinearinterp(x11, x12, x21, x22, y1, y2, c11, c12, c21, c22, x, y, c)

    real(rk), intent(in)  :: x11, x12, x21, x22
    real(rk), intent(in)  :: y1, y2
    real(rk), intent(in)  :: c11, c12, c21, c22
    real(rk), intent(in)  :: x, y
    real(rk), intent(out) :: c
    real(rk)              :: c1, c2

    c1 = (x21 - x) / (x21 - x11) * c11 + (x - x11) / (x21 - x11) * c21
    c2 = (x22 - x) / (x22 - x12) * c12 + (x - x12) / (x22 - x12) * c22
    c = (y2 - y) / (y2 - y1) * c1 + (y - y1) / (y2 - y1) * c2

    return
  end subroutine bilinearinterp
  !===========================================
  subroutine trilinearinterp(x1, x2, y1, y2, z1, z2, &
                             c111, c121, c211, c221, c112, c122, c212, c222, &
                             x, y, z, c)
    !---------------------------------------------
    ! indexing:
    !    c_ijk
    !      111 - left front bottom corner
    !      121 - left back bottom corner
    !      221 - right back bottom corner
    !      222 - right back top corner
    !     c_jk
    !       11 - front bottom corner along x-axis
    !       22 - back top corner along x-axis
    !      c_k
    !        1 - bottom point along y-axis
    !        2 - top point along y-axis
    ! TODO: Maybe calculate dx, dy and dz for each edge?
    !---------------------------------------------

    real(rk), intent(in)  :: x1, x2
    real(rk), intent(in)  :: y1, y2
    real(rk), intent(in)  :: z1, z2
    real(rk)              :: dx, dy, dz
    real(rk), intent(in)  :: c111, c121, c211, c221, c112, c122, c212, c222
    real(rk), intent(in)  :: x, y, z
    real(rk), intent(out) :: c
    real(rk)              :: c11, c12, c21, c22
    real(rk)              :: c1, c2

    dx = (x - x1) / (x2 - x1)
    dy = (y - y1) / (y2 - y1)
    dz = (z - z1) / (z2 - z1)

    c11 = c111 * (1 - dx) + c211 * dx
    c12 = c112 * (1 - dx) + c212 * dx
    c21 = c121 * (1 - dx) + c221 * dx
    c22 = c122 * (1 - dx) + c222 * dx

    c1 = c11 * (1 - dy) + c21 * dy
    c2 = c12 * (1 - dy) + c22 * dy

    c = c1 * (1 - dz) + c2 * dz

    return
  end subroutine trilinearinterp
  !===========================================

end module mod_interp
