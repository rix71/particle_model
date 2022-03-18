#include "cppdefs.h"
module mod_initialise
  !----------------------------------------------------------------
  ! This module is used to read namelists, initialise variables,
  ! allocate arrays etc.
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  use run_params, only: runid, dry_run
  use mod_fieldset
  use field_vars, only: GETMPATH, PMAPFILE, has_subdomains, has_density, has_viscosity, &
                        file_prefix, file_suffix, nlevels, &
                        uvarname, vvarname, wvarname, zaxvarname, elevvarname, zax_style, run_3d, fieldset
  use mod_domain_vars, only: TOPOFILE, bathyvarname, lonvarname, latvarname, nx, ny, domain
  use mod_domain
  ! use fields, only: init_dirlist, init_fields, find_file, find_folder
  use nc_manager, only: nc_read_time_val, nc_var_exists
  use mod_particle_vars, only: inputstep, particle_init_method, coordfile, &
                               max_age, &
                               kill_beached, kill_boundary, &
                               init_particles_from_netcdf, init_particles_from_coordfile
  use time_vars
  use mod_datetime
  use mod_params, only: do_diffusion, do_velocity, Ah, kv
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
    if (ierr .ne. 0) call throw_error('initialise :: init_run', "Failed to open "//NMLFILENAME, ierr)
    read (NMLFILE, nml=run_params)
    close (NMLFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error('initialise :: init_run', "Failed to close "//NMLFILENAME, ierr)

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
      file_prefix, file_suffix, nlevels, &
      uvarname, vvarname, wvarname, zaxvarname, elevvarname, zax_style, run_3d

    FMT1, "======== Init namelist ========"

    open (NMLFILE, file=NMLFILENAME, action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error('initialise :: init_namelist', "Failed to open "//NMLFILENAME, ierr)
    read (NMLFILE, nml=params)
    read (NMLFILE, nml=domain_vars)
    read (NMLFILE, nml=particle_vars)
    read (NMLFILE, nml=time_vars)
    read (NMLFILE, nml=field_vars)
    close (NMLFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error('initialise :: init_namelist', "Failed to close "//NMLFILENAME, ierr)

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
    FMT3, var2val(nlevels)
    FMT3, var2val_char(uvarname)
    FMT3, var2val_char(vvarname)
    FMT3, var2val_char(zaxvarname)
    FMT3, var2val_char(elevvarname)
    FMT3, var2val(zax_style)
    FMT3, var2val(run_3d)

    FMT2, "Finished init namelist"

  end subroutine init_namelist
  !===========================================
  subroutine init_domain
    FMT1, "======== Init domain ========"

    domain = t_domain(nx, ny, topofile=TOPOFILE, lon=lonvarname, lat=latvarname, bathy=bathyvarname)

  end subroutine init_domain
  !===========================================
  subroutine init_time

    character(len=LEN_CHAR_S), parameter :: startName = "Start"
    character(len=LEN_CHAR_S), parameter :: endName = "End"
    character(len=LEN_CHAR_S), parameter :: simName = "Simulation time"
    ! real(rk)                     :: t1, t2

    FMT1, "======== Init time ========"

    ! Initialise datetime (really don't need the names...)
    run_start_dt = t_datetime(run_start)
    call run_start_dt%setName(startName)
    run_end_dt = t_datetime(run_end)
    call run_end_dt%setName(endName)
    theDate = t_datetime(run_start)
    call theDate%setName(simName)

    ! Number of iterations
    nTimes = int(date_diff(run_start_dt, run_end_dt) / dt)

    ! ! Time increment in netCDF input
    ! select case (has_subdomains)
    ! case (.true.)
    !   call find_folder(theDate, initPath)
    !   call nc_read_time_val(trim(initPath)//PROC0, 1, t1)
    !   call nc_read_time_val(trim(initPath)//PROC0, 2, t2)
    ! case (.false.)
    !   call find_file(theDate, initPath)
    !   call nc_read_time_val(trim(initPath), 1, t1)
    !   call nc_read_time_val(trim(initPath), 2, t2)
    ! end select
    ! nc_timestep = t2 - t1

    FMT2, "Model runs from"
    call run_start_dt%print_short_date
    FMT2, "to"
    call run_end_dt%print_short_date
    FMT2, "Timestep: ", dt, " seconds, ", nTimes, " iterations"
    ! FMT2, "netCDF timestep: ", nc_timestep, " seconds"

    FMT2, "Finished init time"

  end subroutine init_time
  !===========================================
  subroutine init_fieldset
    !---------------------------------------------
    ! Allocate arrays for current data.
    ! TODO: Right now it is assumed that all subdomains
    !       are the same size. Also it is assumed that subdomains
    !       exist at all. There should be a switch (e.g. has_subdomains).
    !       This also changes allocation.
    ! EDIT: Only full domain is used, even when the data is in chunks,
    !       since reading into subdomains is not ready yet.
    !---------------------------------------------
    character(len=LEN_CHAR_L) :: initPath
    character(len=LEN_CHAR_L) :: filename
    real(rk)                  :: real_var
    integer(rk)               :: field_count = 0, field_mem

    FMT1, "======== Init fields ========"

    field_mem = (nx * ny * nlevels * 2) * sizeof(real_var)
    FMT2, "Allocating fields of size (nx, ny, nz): (", nx, ", ", ny, ", ", nlevels, ")", &
      field_mem, " bytes per field"
    if (has_subdomains) then
      fieldset = t_fieldset(nx, ny, nlevels, &
                            file_prefix=trim(file_prefix), file_suffix=trim(file_suffix), &
                            domain=domain, &
                            path=GETMPATH, pmap=PMAPFILE)
    else
      fieldset = t_fieldset(nx, ny, nlevels, &
                            file_prefix=trim(file_prefix), file_suffix=trim(file_suffix), &
                            domain=domain, &
                            path=GETMPATH)
    end if

    call fieldset%set_start_time(run_start_dt)

    select case (has_subdomains)
    case (.true.)
      initPath = fieldset%get_folder(1)
      write (filename, '(a)') trim(initPath)//PROC0
    case (.false.)
      initPath = fieldset%get_file(1)
      write (filename, '(a)') trim(initPath)
    end select

    !---------------------------------------------
    ! Fields for h. velocity
    if (nc_var_exists(trim(filename), trim(uvarname))) then
      call fieldset%add_field("U", uvarname)
    else
      call throw_error("initialise :: init_fieldset", "Variable for U velocity does not exist: "//trim(uvarname))
    end if

    if (nc_var_exists(trim(filename), trim(vvarname))) then
      call fieldset%add_field("V", vvarname)
    else
      call throw_error("initialise :: init_fieldset", "Variable for V velocity does not exist: "//trim(vvarname))
    end if
    field_count = field_count + 2
    !---------------------------------------------
    ! Fields for v. velocity
    if (run_3d) then
      if (nc_var_exists(trim(filename), trim(wvarname))) then
        call fieldset%add_field("W", wvarname)
      else
        call throw_error("initialise :: init_fieldset", "Variable for W velocity does not exist: "//trim(wvarname))
      end if

      if (nc_var_exists(trim(filename), trim(zaxvarname))) then
        call fieldset%add_field("ZAX", zaxvarname)
        call fieldset%set_zax("ZAX", zax_style)
      else
        call throw_error("initialise :: init_fieldset", "Variable for Z axis does not exist: "//trim(zaxvarname))
      end if

      field_count = field_count + 2

      if (nc_var_exists(trim(filename), trim(elevvarname))) then
        call fieldset%add_field("ELEV", elevvarname, is_2d=.true.)
        ! set_elev for faster lookup
        field_count = field_count + 1
      end if
    end if
    !---------------------------------------------
    if (do_velocity) then
      !---------------------------------------------
      ! Field for density
      if (nc_var_exists(trim(filename), "rho")) then
        call fieldset%add_field("RHO", "rho")
        has_density = DENSITY
        field_count = field_count + 1
      else
        call throw_warning("initialise :: init_fieldset", "Could not find density ('rho') in "//trim(filename))
        if (nc_var_exists(trim(filename), "temp") .and. &
            nc_var_exists(trim(filename), "salt")) then
          call fieldset%add_field("TEMP", "temp")
          call fieldset%add_field("SALT", "salt")
          has_density = TEMP_SALT
          field_count = field_count + 2
        else
          call throw_warning("initialise :: init_fieldset", "Could not find temperature or salinity ('temp'/'salt') in "//trim(filename)// &
                             ". Using default density.")
          has_density = DEFAULT_DENSITY
        end if
      end if
      !---------------------------------------------
      ! Field for viscosity
      if (nc_var_exists(trim(filename), "nuh")) then
        call fieldset%add_field("VISC", "nuh")
        has_viscosity = .true.
        field_count = field_count + 1
      else
        call throw_warning("initialise :: init_fieldset", "Could not find viscosity ('nuh') in "//trim(filename)// &
                           ". Using default viscosity.")
      end if
    end if
    !---------------------------------------------
    ! if (has_subdomains) call fieldset%init_proc_mask()

    FMT2, "Fields allocated: total", field_count * field_mem, " bytes"

  end subroutine init_fieldset
  !===========================================
  subroutine init_model
    !---------------------------------------------
    ! Call all the subroutines to initialise the model
    !---------------------------------------------

    call init_namelist                 ! init.f90
    call init_domain                   ! init.f90
    call init_time                     ! init.f90
    call init_fieldset                 ! init.f90
    select case (particle_init_method)
    case (TXT_FILE)
      call init_particles_from_coordfile ! particle.f90
    case (NC_FILE)
      call init_particles_from_netcdf    ! particle.f90
    end select

  end subroutine init_model

end module mod_initialise
