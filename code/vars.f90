#include "cppdefs.h"
module run_params
  !----------------------------------------------------------------
  ! ...
  !----------------------------------------------------------------
  use mod_precdefs

  character(len=LEN_CHAR_L) :: runid
  logical                   :: dry_run
  logical                   :: restart
  character(len=LEN_CHAR_L) :: restart_path

end module run_params
!===================================================
module mod_params
  !----------------------------------------------------------------
  ! This includes the model parameters/constants
  ! TODO (later): Diffusion? Biofouling? ...?
  !   viscosity and density defaults in namelist?
  !----------------------------------------------------------------
  use mod_precdefs

  logical             :: do_diffusion, &
                         do_velocity, &                               ! Calculate particles' own velocity
                         do_biofouling, &                             ! Calculate density change due to biological growth
                         run_3d                                       ! 2D or 3D
  integer             :: advection_method
  real(rk)            :: diffusion_hor_const, diffusion_vert_const, & ! Horisontal and vertical diffusion coefs
                         Cm_smagorinsky, &                            ! Used in Ah calculation (Smagorinsky)
                         Im, &                                        ! Light intensity at noon
                         eps_light                                    ! Light extinction coefficient
  real(rk), parameter :: pi = 4.*atan(1.), &                          ! 3, plus a little extra
                         g = 9.81, &
                         k_b = 1.380649e-23, &                        ! Boltzmann constant [m2 kg s-2 K-1]
                         mu_default = 0.0010016, &                    ! Default dynamic viscosity of water at 20 degrees Celcius [N s/m2]
                         rho_default = 1000.0d0                       ! Default density

end module mod_params
!===================================================
module mod_domain_vars
  !----------------------------------------------------------------
  ! This includes variables related to the domain (topography file) (and subdomains?)
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_domain

  integer                   :: nx, ny                ! Size of the domain
  character(len=LEN_CHAR_S) :: bathyvarname, &       ! Bathymetry variable name
                               lonvarname, &         ! Longitude variable name
                               latvarname            ! Latitude variable name
  character(len=LEN_CHAR_L) :: TOPOFILE              ! Topography file
  type(t_domain)            :: domain

end module mod_domain_vars
!===================================================
module time_vars
  !----------------------------------------------------------------
  ! This includes time variables such as:
  ! - model start/end dates
  ! - time step
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_datetime, only: t_datetime

  integer                   :: nTimes                   ! Number of timesteps
  type(t_datetime)          :: theDate, &               ! Current date during simulation
                               run_start_dt, run_end_dt ! Model start/end dates
  character(len=LEN_CHAR_S) :: run_start, run_end       ! Date strings from namelist
  real(rk)                  :: dt, &                    ! Time increment
                               nc_timestep              ! Time increment in input netcdf data

end module time_vars
!===================================================
module field_vars
  !----------------------------------------------------------------
  ! This includes variables related to the hydrodynamic data
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_fieldset

  logical                   :: has_subdomains, &               ! Is the data in multiple files (true) or one file (false)?
                               has_viscosity
  integer                   :: nlevels, &
                               zax_style, &                    ! Depth values (1) or layer thickness (2)
                               has_density                     ! 0 - default density, 1 - has variable, 2 - calculate from T/S
  character(len=LEN_CHAR_S) :: uvarname, vvarname, wvarname, & ! Names of the variables. Necessary?
                               zaxvarname, elevvarname, &
                               rhovarname, tempvarname, &
                               saltvarname, viscvarname
  character(len=LEN_CHAR_L) :: file_prefix, file_suffix        ! What comes before and after the proc. number?
  character(len=LEN_CHAR_L) :: GETMPATH, PMAPFILE              ! Path to GETM output and processor map
  type(t_fieldset)          :: fieldset

end module field_vars
