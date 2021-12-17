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

  logical             :: do_diffusion
  logical             :: do_velocity      ! Calculate particles' own velocity
  real(rk)            :: Ah, kv           ! Horisontal and vertical diffusion coefs
  real(rk), parameter :: pi = 4.*atan(1.) ! 3, plus a little extra
  real(rk), parameter :: g = 9.81
  real(rk), parameter :: mu = 0.0010016   ! Dynamic viscosity of water at 20 degrees Celcius [N s/m2]

end module params
!===================================================
module loop_vars
  !----------------------------------------------------------------
  ! This module contains temporary variables used in the loop
  ! ??? Is this even necessary? ???
  ! As it turns out, yes.
  !----------------------------------------------------------------
  use precdefs
  use modtime, only: datetime

  integer            :: pnum                 ! Processor/subdomain number
  integer            :: ig, jg, kg           ! Global indices
  integer            :: il, jl, kl           ! Local indices (subdomain) TODO: Probably won't need local k indices...
  integer            :: ncNTimes             ! Number of timesteps in nc file
  integer            :: nc_itime             ! Time index in netCDF
  integer            :: nc_itime_next        ! Next time index in netCDF
  type(datetime)     :: dateThis, dateNext   ! Temporary date variables, first date of netCDF files
  type(datetime)     :: dateThisNcTimestep   ! Date of current timestep in netCDF
  type(datetime)     :: dateNextNcTimestep   ! Date of next timestep in netCDF
  character(len=256) :: thisPath, nextPath   ! Temporary path to data
  real(rk)           :: ilr, jlr, klr        ! Local real indices for interpolation
  real(rk)           :: igr, jgr, kgr        ! Global real indices for interpolation
  real(rk)           :: xnow, xnew           ! Temporary location in lon, lat
  real(rk)           :: ynow, ynew
  real(rk)           :: znow, znew           ! Vertical position
  real(rk)           :: cartx, cartx_new     ! Temporary loc. in Cart. coordinates
  real(rk)           :: carty, carty_new
  real(rk)           :: pvelu, pvelunew      ! Temporary particle velocity
  real(rk)           :: pvelv, pvelvnew      ! Temporary particle velocity
  real(rk)           :: pvelw, pvelwnew      ! Temporary particle velocity

end module loop_vars
!===================================================
module domain_vars
  !----------------------------------------------------------------
  ! This includes variables related to the domain (topography file) (and subdomains?)
  !----------------------------------------------------------------
  use precdefs

  integer                :: nx, ny             ! Size of the domain
  integer, allocatable   :: seamask(:, :)       ! Global seamask
  character(len=128)     :: bathyvarname       ! Bathymetry variable name
  character(len=128)     :: lonvarname         ! Longitude variable name
  character(len=128)     :: latvarname         ! Latitude variable name
  character(len=512)     :: TOPOFILE           ! Topography file
  real(rk)               :: x0, x1, y0, y1     ! Minimum and maximum coordinates
  real(rk)               :: dx, dy, dx_m, dy_m ! Spatial steps (lon/lat and meters)
  real(rk)               :: dz                 ! Vertical step (meters)
  real(rk), allocatable  :: lons(:), lats(:)   ! Global longitude/latitude
  real(rk), allocatable  :: depdata(:, :)       ! Global bathymetry

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
  type(datetime)    :: theDate                  ! Current date during simulation
  type(datetime)    :: run_start_dt, run_end_dt ! Model start/end dates
  character(len=64) :: run_start, run_end       ! Date strings from namelist
  real(rk)          :: dt                       ! Time increment
  real(rk)          :: nc_timestep              ! Time increment in input netcdf data

end module time_vars
!===================================================
module field_vars
  !----------------------------------------------------------------
  ! This includes variables related to the (current) field data (H.D. model output file(s))
  !----------------------------------------------------------------
  use precdefs

  logical               :: has_subdomains                    ! Is the data in multiple files (true) or one file (false)?
  logical               :: run_3d                            ! 2D or 3D
  integer               :: startlevel, nlevels
  integer               :: zax_style                         ! Depth values (1) or layer thickness (2)
  character(len=32)     :: uvarname, vvarname, wvarname      ! Names of the variables. Necessary?
  character(len=32)     :: zaxvarname                        ! Names of the variables. Necessary?
  character(len=128)    :: file_prefix, file_suffix          ! What comes before and after the proc. number?
  character(len=512)    :: GETMPATH, PMAPFILE                ! Path to GETM output and processor map
  real(rk)              :: rho_sea                           ! Seawater density
  real(rk)              :: viscosity                         ! Viscosity
  real(rk)              :: uspeed, vspeed, wspeed            ! Speed at the location of the particle
  real(rk)              :: uspeednew, vspeednew, wspeednew
  real(rk), allocatable :: udata(:, :, :), vdata(:, :, :)        ! The full current field of the (sub)domain
  real(rk), allocatable :: wdata(:, :, :)                        ! The full current field of the (sub)domain
  real(rk), allocatable :: udatanew(:, :, :), vdatanew(:, :, :)  ! Same thing, but at the next timestep
  real(rk), allocatable :: wdatanew(:, :, :)                     ! Same thing, but at the next timestep
  real(rk), allocatable :: zaxdata(:, :, :), zaxdatanew(:, :, :) ! Layer thickness or depth (TODO: should distinguish the two!!!)
  real(rk), allocatable :: udata_interp(:, :, :), vdata_interp(:, :, :), wdata_interp(:, :, :)
  real(rk), allocatable :: udatanew_interp(:, :, :), vdatanew_interp(:, :, :), wdatanew_interp(:, :, :)
  real(rk), allocatable :: zaxdata_interp(:, :, :), zaxdatanew_interp(:, :, :)

end module field_vars
