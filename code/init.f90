#include "cppdefs.h"
module mod_initialise
  !----------------------------------------------------------------
  ! This module is used to read namelists, initialise variables,
  ! allocate arrays etc.
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  use run_params, only: runid, dry_run, restart, restart_path
  use mod_fieldset
  use field_vars, only: GETMPATH, PMAPFILE, has_subdomains, has_density, has_viscosity, &
                        file_prefix, file_suffix, nlevels, &
                        uvarname, vvarname, wvarname, zaxvarname, elevvarname, rhovarname, &
                        tempvarname, saltvarname, viscvarname, zax_style, fieldset
  use mod_domain_vars, only: TOPOFILE, bathyvarname, lonvarname, latvarname, nx, ny, domain
  use mod_domain
  use nc_manager, only: nc_read_time_val, nc_var_exists
  use mod_particle_vars, only: inputstep, particle_init_method, coordfile, &
                               max_age, &
                               kill_beached, kill_boundary, &
                               init_particles
  use time_vars
  use mod_datetime
  use mod_biofouling, only: init_biofouling
  use mod_params, only: do_diffusion, do_velocity, do_biofouling, advection_method, &
                        diffusion_hor_const, diffusion_vert_const, run_3d, Cm_smagorinsky
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
    namelist /run_params/ runid, dry_run, restart, restart_path

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
    namelist /params/ do_diffusion, do_velocity, do_biofouling, run_3d, advection_method, diffusion_hor_const, diffusion_vert_const, Cm_smagorinsky
    namelist /domain_vars/ TOPOFILE, bathyvarname, lonvarname, latvarname, nx, ny
    namelist /particle_vars/ inputstep, particle_init_method, coordfile, max_age, kill_beached, kill_boundary
    namelist /time_vars/ run_start, run_end, dt
    namelist /field_vars/ GETMPATH, PMAPFILE, has_subdomains, &
      file_prefix, file_suffix, nlevels, &
      uvarname, vvarname, wvarname, zaxvarname, elevvarname, rhovarname, &
      tempvarname, saltvarname, viscvarname, zax_style

    FMT1, "======== Init namelist ========"

    open (NMLFILE, file=NMLFILENAME, action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error('initialise :: init_namelist', "Failed to open "//NMLFILENAME, ierr)
    read (NMLFILE, nml=params, iostat=ierr); if (ierr .ne. 0) call throw_error("initialise :: init_namelist", "Could not read 'params' namelist", ierr)
    read (NMLFILE, nml=domain_vars, iostat=ierr); if (ierr .ne. 0) call throw_error("initialise :: init_namelist", "Could not read 'domain_vars' namelist", ierr)
    read (NMLFILE, nml=particle_vars, iostat=ierr); if (ierr .ne. 0) call throw_error("initialise :: init_namelist", "Could not read 'particle_vars' namelist", ierr)
    read (NMLFILE, nml=time_vars, iostat=ierr); if (ierr .ne. 0) call throw_error("initialise :: init_namelist", "Could not read 'time_vars' namelist", ierr)
    read (NMLFILE, nml=field_vars, iostat=ierr); if (ierr .ne. 0) call throw_error("initialise :: init_namelist", "Could not read 'field_vars' namelist", ierr)
    close (NMLFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error('initialise :: init_namelist', "Failed to close "//NMLFILENAME, ierr)

    FMT2, LINE
    FMT2, "&params"
    FMT3, var2val(do_diffusion)
    FMT3, var2val(do_velocity)
    FMT3, var2val(do_biofouling)
    FMT3, var2val(run_3d)
    FMT3, var2val(advection_method)
    FMT3, var2val(diffusion_hor_const)
    FMT3, var2val(diffusion_vert_const)
    FMT3, var2val(Cm_smagorinsky)
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
    FMT3, var2val_char(rhovarname)
    FMT3, var2val_char(tempvarname)
    FMT3, var2val_char(saltvarname)
    FMT3, var2val_char(viscvarname)
    FMT3, var2val(zax_style)

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

    FMT2, "Model runs from"
    call run_start_dt%print_short_date
    FMT2, "to"
    call run_end_dt%print_short_date
    FMT2, "Timestep: ", dt, " seconds, ", nTimes, " iterations"

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
    FMT2, "Using full domain"
    FMT2, "Allocating fields of size (nx, ny, nz): (", nx, ", ", ny, ", ", nlevels, ")", field_mem, " bytes per field"
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
    call fieldset%set_simulation_timestep(dt)

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
      call fieldset%set_u_component("U")
    else
      call throw_error("initialise :: init_fieldset", "Variable for U velocity does not exist: "//trim(uvarname))
    end if

    if (nc_var_exists(trim(filename), trim(vvarname))) then
      call fieldset%add_field("V", vvarname)
      call fieldset%set_v_component("V")
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
      if (nc_var_exists(trim(filename), trim(rhovarname))) then
        call fieldset%add_field("RHO", rhovarname)
        has_density = DENSITY
        field_count = field_count + 1
      else
        call throw_warning("initialise :: init_fieldset", "Could not find density ('rho') in "//trim(filename))
        if (nc_var_exists(trim(filename), trim(tempvarname)) .and. &
            nc_var_exists(trim(filename), trim(saltvarname))) then
          call fieldset%add_field("TEMP", tempvarname)
          call fieldset%add_field("SALT", saltvarname)
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
      if (nc_var_exists(trim(filename), trim(viscvarname))) then
        call fieldset%add_field("VISC", viscvarname)
        has_viscosity = .true.
        field_count = field_count + 1
      else
        call throw_warning("initialise :: init_fieldset", "Could not find viscosity ('nuh') in "//trim(filename)// &
                           ". Using default viscosity.")
      end if
    end if

    !---------------------------------------------
    if (do_biofouling) then
      call init_biofouling(fieldset, filename)
    end if
    !---------------------------------------------
    ! if (has_subdomains) call fieldset%init_proc_mask()

    FMT1, LINE
    FMT2, "Fields allocated: total", fieldset%num_fields * field_mem, " bytes"
    call fieldset%list_fields()

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
    call init_particles

  end subroutine init_model

end module mod_initialise
