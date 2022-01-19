#include "cppdefs.h"
module run_params
  !----------------------------------------------------------------
  ! ...
  !----------------------------------------------------------------

  character(len=128) :: runid
  logical            :: dry_run

end module run_params
!===================================================
module params
  !----------------------------------------------------------------
  ! This includes the model parameters/constants
  ! TODO (later): Diffusion? Biofouling? ...?
  !----------------------------------------------------------------
  use precdefs

  logical             :: do_diffusion, &
                         do_velocity         ! Calculate particles' own velocity
  real(rk)            :: Ah, kv              ! Horisontal and vertical diffusion coefs
  real(rk), parameter :: pi = 4.*atan(1.), & ! 3, plus a little extra
                         g = 9.81, &
                         mu = 0.0010016      ! Dynamic viscosity of water at 20 degrees Celcius [N s/m2]

end module params
!===================================================
module loop_vars
  !----------------------------------------------------------------
  ! This module contains temporary variables used in the loop
  !----------------------------------------------------------------
  use precdefs
  use modtime, only: datetime

  integer            :: pnum, &               ! Processor/subdomain number
                        ig, jg, kg, &         ! Global indices
                        ncNTimes, &           ! Number of timesteps in nc file
                        nc_itime, &           ! Time index in netCDF
                        nc_itime_next         ! Next time index in netCDF
  type(datetime)     :: dateThis, dateNext, & ! Temporary date variables, first date of netCDF files
                        dateThisNcTimestep, & ! Date of current timestep in netCDF
                        dateNextNcTimestep    ! Date of next timestep in netCDF
  character(len=256) :: thisPath, nextPath    ! Temporary path to data
  real(rk)           :: igr, jgr, kgr, &      ! Global real indices for interpolation
                        xnow, xnew, &         ! Temporary location in lon, lat
                        ynow, ynew, &
                        znow, znew, &         ! Vertical position
                        cartx, cartx_new, &   ! Temporary loc. in Cart. coordinates
                        carty, carty_new, &
                        pvelu, pvelunew, &    ! Temporary particle velocity
                        pvelv, pvelvnew, &    ! Temporary particle velocity
                        pvelw, pvelwnew       ! Temporary particle velocity

end module loop_vars
!===================================================
module domain_vars
  !----------------------------------------------------------------
  ! This includes variables related to the domain (topography file) (and subdomains?)
  !----------------------------------------------------------------
  use precdefs

  integer                :: nx, ny                ! Size of the domain
  integer, allocatable   :: seamask(:, :)         ! Global seamask
  character(len=128)     :: bathyvarname, &       ! Bathymetry variable name
                            lonvarname, &         ! Longitude variable name
                            latvarname            ! Latitude variable name
  character(len=512)     :: TOPOFILE              ! Topography file
  real(rk)               :: x0, x1, y0, y1, &     ! Minimum and maximum coordinates
                            dx, dy, dx_m, dy_m, & ! Spatial steps (lon/lat and meters)
                            dz                    ! Vertical step (meters)
  real(rk), allocatable  :: lons(:), lats(:), &   ! Global longitude/latitude
                            depdata(:, :)         ! Global bathymetry

end module domain_vars
!===================================================
module time_vars
  !----------------------------------------------------------------
  ! This includes time variables such as:
  ! - model start/end dates
  ! - time step
  !----------------------------------------------------------------
  use precdefs
  use modtime, only: datetime

  integer           :: nTimes                   ! Number of timesteps
  type(datetime)    :: theDate, &               ! Current date during simulation
                       run_start_dt, run_end_dt ! Model start/end dates
  character(len=64) :: run_start, run_end       ! Date strings from namelist
  real(rk)          :: dt, &                    ! Time increment
                       nc_timestep              ! Time increment in input netcdf data

end module time_vars
!===================================================
module field_vars
  !----------------------------------------------------------------
  ! This includes variables related to the hydrodynamic data
  !----------------------------------------------------------------
  use precdefs

  logical                                   :: has_subdomains, &               ! Is the data in multiple files (true) or one file (false)?
                                               has_viscosity, &
                                               run_3d                          ! 2D or 3D
  integer                                   :: startlevel, nlevels, &
                                               zax_style, &                    ! Depth values (1) or layer thickness (2)
                                               has_density                     ! 0 - default density, 1 - has variable, 2 - calculate from T/S
  character(len=32)                         :: uvarname, vvarname, wvarname, & ! Names of the variables. Necessary?
                                               zaxvarname                      ! Names of the variables. Necessary?
  character(len=128)                        :: file_prefix, file_suffix        ! What comes before and after the proc. number?
  character(len=512)                        :: GETMPATH, PMAPFILE              ! Path to GETM output and processor map
  real(rk)                                  :: rho_sea, &                      ! Seawater density
                                               viscosity, &                    ! Viscosity
                                               uspeed, vspeed, wspeed, &       ! Speed at the location of the particle
                                               uspeednew, vspeednew, wspeednew
  real(rk), dimension(:, :, :), allocatable, target :: udata, vdata, &                 ! The full current field of the (sub)domain
                                                       wdata, &                        ! The full current field of the (sub)domain
                                                       udatanew, vdatanew, &           ! Same thing, but at the next timestep
                                                       wdatanew, &                     ! Same thing, but at the next timestep
                                                       zaxdata, zaxdatanew, &          ! Layer thickness or depth (TODO: should distinguish the two!!!)
                                                       udata_interp, vdata_interp, wdata_interp, &
                                                       udatanew_interp, vdatanew_interp, wdatanew_interp, &
                                                       zaxdata_interp, zaxdatanew_interp, &
                                                       visc, viscnew, &                ! Viscosity
                                                       density, densitynew, &          ! Density
                                                       temp, tempnew, &                ! Temperature
                                                       salt, saltnew                   ! Salinity

end module field_vars
