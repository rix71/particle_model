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
#ifdef DEBUG
    use nc_manager
    integer :: nc_x_dimid, nc_y_dimid
#endif

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

    FMT2, LINE
    FMT2, "Coordinates:"
    FMT3, var2val(y0), "[deg N], ", var2val(y1), "[deg N]"
    FMT3, var2val(x0), "[deg E], ", var2val(x1), "[deg E]"
    FMT2, "Cell size:"
    FMT3, var2val(dy), "[deg], ", var2val(dy_m), "[m]"
    FMT3, var2val(dx), "[deg], ", var2val(dx_m), "[m]"

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
    do ii = 1, nx
      if (depdata(ii, 1) .gt. 0.0d0) then
        seamask(ii, 1) = 4
      else
        seamask(ii, 1) = 1
      end if
      if (depdata(ii, ny) .gt. 0.0d0) then
        seamask(ii, ny) = 4
      else
        seamask(ii, ny) = 1
      end if
    end do
    do jj = 1, ny
      if (depdata(1, jj) .gt. 0.0d0) then
        seamask(1, jj) = 4
      else
        seamask(1, jj) = 1
      end if
      if (depdata(nx, jj) .gt. 0.0d0) then
        seamask(nx, jj) = 4
      else
        seamask(nx, jj) = 1
      end if
    end do

    FMT2, sum(seamask, mask=seamask == 1), " land points"
    FMT2, sum(seamask, mask=seamask == 2) / 2, " sea points"
    FMT2, sum(seamask, mask=seamask == 3) / 3, " beach points"
    FMT2, sum(seamask, mask=seamask == 4) / 4, " boundary points"
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
    call nc_write(FNAME, seamask, "seamask", nx, ny)
    call nc_write(FNAME, lons, "lon", nx)
    call nc_write(FNAME, lats, "lat", ny)
#undef FNAME
#endif

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
