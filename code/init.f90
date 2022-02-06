#include "cppdefs.h"
module initialise
  !----------------------------------------------------------------
  ! This module is used to read namelists, initialise variables,
  ! allocate arrays etc.
  !----------------------------------------------------------------
  use precdefs
  use errors
  use run_params, only: runid, dry_run
  use field_vars, only: GETMPATH, PMAPFILE, has_subdomains, &
                        file_prefix, file_suffix, startlevel, nlevels, &
                        uvarname, vvarname, wvarname, zaxvarname, zax_style, run_3d
  use domain_vars, only: TOPOFILE, bathyvarname, lonvarname, latvarname, nx, ny, lons, lats
  use domain, only: init_domain
  use fields, only: init_dirlist, init_fields, find_file, find_folder
  use nc_manager, only: nc_read_time_val
  use particle_vars, only: inputstep, particle_init_method, coordfile, &
                           max_age, &
                           kill_beached, kill_boundary, &
                           init_particles_from_netcdf, init_particles_from_coordfile
  use time_vars
  use modtime
  use params, only: do_diffusion, do_velocity, Ah, kv
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: init_run, init_model
  !---------------------------------------------
  integer :: ierr
  !===================================================
contains
  !===========================================
  subroutine init_run
    !---------------------------------------------
    ! First things to read in
    !---------------------------------------------
    namelist /run_params/ runid, dry_run

    open (NMLFILE, file=NMLFILENAME, action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error('init_run', "Failed to open "//NMLFILENAME, ierr)
    read (NMLFILE, nml=run_params)
    close (NMLFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error('init_run', "Failed to close "//NMLFILENAME, ierr)

  end subroutine init_run
  !===========================================
  subroutine init_namelist
    !---------------------------------------------
    ! Namelists
    !---------------------------------------------
    namelist /params/ do_diffusion, do_velocity, Ah, kv
    namelist /domain_vars/ TOPOFILE, bathyvarname, lonvarname, latvarname, nx, ny
    namelist /particle_vars/ inputstep, particle_init_method, coordfile, max_age, kill_beached, kill_boundary
    namelist /time_vars/ run_start, run_end, dt
    namelist /field_vars/ GETMPATH, PMAPFILE, has_subdomains, &
      file_prefix, file_suffix, startlevel, nlevels, &
      uvarname, vvarname, wvarname, zaxvarname, zax_style, run_3d

    FMT1, "======== Init namelist ========"

    open (NMLFILE, file=NMLFILENAME, action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error('init_namelist', "Failed to open "//NMLFILENAME, ierr)
    read (NMLFILE, nml=params)
    read (NMLFILE, nml=domain_vars)
    read (NMLFILE, nml=particle_vars)
    read (NMLFILE, nml=time_vars)
    read (NMLFILE, nml=field_vars)
    close (NMLFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error('init_namelist', "Failed to close "//NMLFILENAME, ierr)

    FMT2, LINE
    FMT2, "&params"
    FMT3, var2val(do_diffusion)
    FMT3, var2val(do_velocity)
    FMT3, var2val(Ah)
    FMT3, var2val(kv)
    FMT2, LINE
    FMT2, "&domain_vars"
    FMT3, var2val_char(TOPOFILE)
    FMT3, var2val_char(bathyvarname)
    FMT3, var2val_char(lonvarname)
    FMT3, var2val_char(latvarname)
    FMT3, var2val(nx)
    FMT3, var2val(ny)
    FMT2, LINE
    FMT2, "&particle_vars"
    FMT3, var2val(inputstep)
    FMT3, var2val(particle_init_method)
    FMT3, var2val_char(coordfile)
    FMT3, var2val(max_age)
    FMT3, var2val(kill_beached)
    FMT3, var2val(kill_boundary)
    FMT2, LINE
    FMT2, "&field_vars"
    FMT3, var2val_char(GETMPATH)
    FMT3, var2val_char(PMAPFILE)
    FMT3, var2val(has_subdomains)
    FMT3, var2val_char(file_prefix)
    FMT3, var2val_char(file_suffix)
    FMT3, var2val(startlevel)
    FMT3, var2val(nlevels)
    FMT3, var2val_char(uvarname)
    FMT3, var2val_char(vvarname)
    FMT3, var2val_char(zaxvarname)
    FMT3, var2val(zax_style)
    FMT3, var2val(run_3d)

    FMT2, "Finished init namelist"

  end subroutine init_namelist
  !===========================================
  subroutine init_time

    character(len=16), parameter :: startName = "Start"
    character(len=16), parameter :: endName = "End"
    character(len=16), parameter :: simName = "Simulation time"
    character(len=256)           :: initPath
    real(rk)                     :: t1, t2

    FMT1, "======== Init time ========"

    ! Initialise datetime (really don't need the names...)
    run_start_dt = t_datetime(run_start)
    call run_start_dt%setName(startName)
    run_end_dt = t_datetime(run_end)
    call run_end_dt%setName(endName)
    theDate = t_datetime(run_start)
    call theDate%setName(simName)

    ! Number of iterations
    nTimes = int(dateDiff(run_start_dt, run_end_dt) / dt)

    ! Time increment in netCDF input
    select case (has_subdomains)
    case (.true.)
      call find_folder(theDate, initPath)
      call nc_read_time_val(trim(initPath)//PROC0, 1, t1)
      call nc_read_time_val(trim(initPath)//PROC0, 2, t2)
    case (.false.)
      call find_file(theDate, initPath)
      call nc_read_time_val(trim(initPath), 1, t1)
      call nc_read_time_val(trim(initPath), 2, t2)
    end select
    nc_timestep = t2 - t1

    FMT2, "Model runs from"
    call run_start_dt%print_short_date
    FMT2, "to"
    call run_end_dt%print_short_date
    FMT2, "Timestep: ", dt, " seconds, ", nTimes, " iterations"
    FMT2, "netCDF timestep: ", nc_timestep, " seconds"

    FMT2, "Finished init time"

  end subroutine init_time
  !===========================================
  subroutine init_model
    !---------------------------------------------
    ! Call all the subroutines to initialise the model
    !---------------------------------------------

    call init_namelist                 ! init.f90
    call init_dirlist                  ! fields.f90
    call init_time                     ! init.f90
    call init_domain                   ! domain.f90
    select case (particle_init_method)
    case (TXT_FILE)
      call init_particles_from_coordfile ! particle.f90
    case (NC_FILE)
      call init_particles_from_netcdf    ! particle.f90
    end select
    call init_fields                   ! fields.f90

  end subroutine init_model

end module initialise
