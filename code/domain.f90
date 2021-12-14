#include "cppdefs.h"
module domain
  !----------------------------------------------------------------
  ! Initialise domain and seamask
  !----------------------------------------------------------------
  use precdefs
  use errors
  use params, only: pi
  use domain_vars
  use nc_manager, only: nc_read1d, nc_read2d
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: init_domain, lonlat2xy, xy2lonlat
  !===================================================
contains
  !===========================================
  subroutine init_domain
    !---------------------------------------------
    ! Initialize the global longitude/latitude
    ! and seamask
    !---------------------------------------------

    integer :: ii, jj

    dbghead(init_domain)

    !print
    FMT1, "======== Init domain ========"

    allocate (lons(nx), lats(ny))

    FMT2, "Reading coordinates"
    call nc_read1d(trim(TOPOFILE), trim(lonvarname), nx, lons)
    call nc_read1d(trim(TOPOFILE), trim(latvarname), ny, lats)

    y0 = lats(1); y1 = lats(ny); dy = lats(2) - lats(1); dy_m = dy * 60.*1852.
    x0 = lons(1); x1 = lons(nx); dx = lons(2) - lons(1); dx_m = dx * 60.*1852.*cos(0.5 * (y0 + y1) * pi / 180.)
    debug(y0); debug(y1); debug(dy); 
    debug(dy_m)
    debug(x0); debug(x1); debug(dx); 
    debug(dx_m)

    allocate (depdata(nx, ny), seamask(nx, ny))

    FMT2, "Reading bathymetry"
    call nc_read2d(trim(TOPOFILE), trim(bathyvarname), nx, ny, depdata)

    !---------------------------------------------
    ! TODO: Seamask could have another value (4) to represent boundaries.
    !       Boundary should have a thickness!
    FMT2, "Making seamask"
    do ii = 2, nx - 1
      do jj = 2, ny - 1
        if (depdata(ii, jj) .gt. 0.0d0) then
          if ((depdata(ii + 1, jj) .le. 0.0) .or. (depdata(ii - 1, jj) .le. 0.0) .or. &
              (depdata(ii, jj + 1) .le. 0.0) .or. (depdata(ii, jj - 1) .le. 0.0) .or. &
              (depdata(ii + 1, jj + 1) .le. 0.0) .or. (depdata(ii + 1, jj - 1) .le. 0.0) .or. &
              (depdata(ii - 1, jj - 1) .le. 0.0) .or. (depdata(ii - 1, jj + 1) .le. 0.0)) then
            seamask(ii, jj) = 3
          else
            seamask(ii, jj) = 2
          end if
        else
          seamask(ii, jj) = 1
        end if
      end do
    end do
    seamask(1, :) = 1
    seamask(nx, :) = 1
    seamask(:, 1) = 1

    FMT2, "Finished init_domain"

    dbgtail(init_domain)
    return
  end subroutine init_domain
  !===========================================
  subroutine lonlat2xy(lonin, latin, x_ref, y_ref, xout, yout)

    real(rk), intent(in)  :: lonin, latin
    real(rk), intent(in)  :: x_ref, y_ref
    real(rk), intent(out) :: xout, yout

    xout = (lonin - x_ref) / dx * dx_m
    yout = (latin - y_ref) / dy * dy_m

  end subroutine lonlat2xy
  !===========================================
  subroutine xy2lonlat(xin, yin, x_ref, y_ref, lonout, latout)

    real(rk), intent(in)  :: xin, yin
    real(rk), intent(in)  :: x_ref, y_ref
    real(rk), intent(out) :: lonout, latout

    lonout = xin / dx_m * dx + x_ref
    latout = yin / dy_m * dy + y_ref

  end subroutine xy2lonlat
end module domain
