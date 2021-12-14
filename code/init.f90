#include "cppdefs.h"
module initialise
  !----------------------------------------------------------------
  ! This module is used to read namelists, initialise variables,
  ! allocate arrays etc.
  !----------------------------------------------------------------
  use precdefs
  use errors
  use run_params, only: runid
  use field_vars, only: GETMPATH, PMAPFILE, has_subdomains, &
                        file_prefix, file_suffix, startlevel, nlevels, &
                        uvarname, vvarname, wvarname, zaxvarname, zax_style, run_3d
  use domain_vars, only: TOPOFILE, bathyvarname, lonvarname, latvarname, nx, ny, lons, lats
  use domain, only: init_domain
  use fields, only: init_dirlist, init_fields, find_file, find_folder
  use nc_manager, only: nc_read_time_val
  use particle_vars, only: inputstep, coordfile, nParticles, &
                           nInitParticles, initCoords, particles
  use time_vars
  use modtime
  use loop_vars, only: ipart
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
    namelist /run_params/ runid

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
    ! TODO: Is run_3d in the right place? (Run options, maybe?)
    !---------------------------------------------
    namelist /params/ do_diffusion, do_velocity, Ah, kv
    namelist /domain_vars/ TOPOFILE, bathyvarname, lonvarname, latvarname, nx, ny
    namelist /particle_vars/ inputstep, coordfile
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
    run_start_dt = datetime(run_start)
    call run_start_dt%setName(startName)
    run_end_dt = datetime(run_end)
    call run_end_dt%setName(endName)
    theDate = datetime(run_start)
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
  subroutine init_particles_from_coordfile
    !---------------------------------------------
    ! Allocate array for estimated amount of particles
    !---------------------------------------------

    !print
    FMT1, "======== Init particles ========"

    open (COORDFILE, file=trim(coordfile), action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error("init_particles_from_coordfile", "Failed to open "//trim(coordfile), ierr)
    read (COORDFILE, *) nInitParticles
    allocate (initCoords(nInitParticles, 6))
    do ipart = 1, nInitParticles
      read (COORDFILE, *) initCoords(ipart, 1), initCoords(ipart, 2), &
        initCoords(ipart, 3), initCoords(ipart, 4), &
        initCoords(ipart, 5), initCoords(ipart, 6)
    end do
    close (COORDFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error("init_particles_from_coordfile", "Failed to close "//trim(coordfile), ierr)

    do ipart = 1, nInitParticles
      if (initCoords(ipart, 1) < lons(1) .or. initCoords(ipart, 1) > lons(nx) .or. &
          initCoords(ipart, 2) < lats(1) .or. initCoords(ipart, 2) > lats(ny)) then
        ERROR, "Particle", ipart, ":", initCoords(ipart, 1), initCoords(ipart, 2)
        call throw_error("init_particles_from_coordfile", "Particle initialized outside of domain")
      end if
    end do
    FMT2, "Initial coordinates OK"

    nParticles = nInitParticles * (nTimes / inputstep) + nInitParticles
    allocate (particles(nParticles))
    FMT2, "Allocated array for", nParticles, "particles"

    !print
    FMT2, "Finished init particles"

  end subroutine init_particles_from_coordfile
  !===========================================
  subroutine init_model
    !---------------------------------------------
    ! Call all the subroutines to initialise the model
    ! TODO: options, for example, different ways
    !       to init particles
    !---------------------------------------------

    call init_namelist                 ! init.f90
    call init_dirlist                  ! fields.f90
    call init_time                     ! init.f90
    call init_domain                   ! domain.f90
    call init_particles_from_coordfile ! init.f90
    call init_fields                   ! fields.f90

  end subroutine init_model

end module initialise
