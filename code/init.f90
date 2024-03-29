#include "cppdefs.h"
#include "field.h"
#include "file.h"
module mod_initialise
  !----------------------------------------------------------------
  ! This module is used to read namelists, initialise variables,
  ! allocate arrays etc.
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  use run_params, only: runid, dry_run, restart, restart_path, nmlfilename
  use mod_fieldset
  use field_vars, only: GETMPATH, PMAPFILE, has_subdomains, density_method, viscosity_method, has_bottom_stress, &
                        file_prefix, file_suffix, nlevels, &
                        xdimname, ydimname, zdimname, &
                        uvarname, vvarname, wvarname, zaxvarname, elevvarname, rhovarname, &
                        tempvarname, saltvarname, viscvarname, taubxvarname, taubyvarname, zax_style, zax_direction, fieldset
  use mod_domain_vars, only: TOPOFILE, bathyvarname, lonvarname, latvarname, nx, ny, domain
  use mod_domain
  use nc_manager, only: nc_read_time_val, nc_var_exists
  use mod_particle_vars, only: inputstep, particle_init_method, coordfile, &
                               max_age, kill_beached, kill_boundary, &
                               initialise_particles
  use time_vars
  use mod_datetime
  use mod_biofouling, only: init_biofouling
  use mod_params, only: do_diffusion, do_velocity, do_biofouling, advection_method, &
     diffusion_hor_const, diffusion_vert_const, run_3d, Cm_smagorinsky, resuspension_coeff, resuspension_threshold, roughness_height
  use mod_output, only: init_output
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

    open (NMLFILE, file=trim(nmlfilename), action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error('initialise :: init_run', "Failed to open "//trim(nmlfilename), ierr)
    read (NMLFILE, nml=run_params)
    close (NMLFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error('initialise :: init_run', "Failed to close "//trim(nmlfilename), ierr)

  end subroutine init_run
  !===========================================
  subroutine init_namelist()
    !---------------------------------------------
    ! Namelists
    !---------------------------------------------
    namelist /params/ do_diffusion, do_velocity, do_biofouling, run_3d, advection_method, diffusion_hor_const, diffusion_vert_const, Cm_smagorinsky, resuspension_coeff, resuspension_threshold, roughness_height
    namelist /domain_vars/ TOPOFILE, bathyvarname, lonvarname, latvarname, nx, ny
    namelist /particle_vars/ inputstep, particle_init_method, coordfile, max_age, kill_beached, kill_boundary
    namelist /time_vars/ run_start, run_end, dt
    namelist /field_vars/ GETMPATH, PMAPFILE, has_subdomains, &
      file_prefix, file_suffix, nlevels, &
      xdimname, ydimname, zdimname, &
      uvarname, vvarname, wvarname, zaxvarname, elevvarname, rhovarname, &
      tempvarname, saltvarname, viscvarname, taubxvarname, taubyvarname, zax_style, zax_direction

    FMT1, "======== Init namelist ========"

    open (NMLFILE, file=trim(nmlfilename), action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error('initialise :: init_namelist', "Failed to open "//trim(nmlfilename), ierr)
    read (NMLFILE, nml=params, iostat=ierr); if (ierr .ne. 0) call throw_error("initialise :: init_namelist", "Could not read 'params' namelist", ierr)
    read (NMLFILE, nml=domain_vars, iostat=ierr); if (ierr .ne. 0) call throw_error("initialise :: init_namelist", "Could not read 'domain_vars' namelist", ierr)
    read (NMLFILE, nml=particle_vars, iostat=ierr); if (ierr .ne. 0) call throw_error("initialise :: init_namelist", "Could not read 'particle_vars' namelist", ierr)
    read (NMLFILE, nml=time_vars, iostat=ierr); if (ierr .ne. 0) call throw_error("initialise :: init_namelist", "Could not read 'time_vars' namelist", ierr)
    read (NMLFILE, nml=field_vars, iostat=ierr); if (ierr .ne. 0) call throw_error("initialise :: init_namelist", "Could not read 'field_vars' namelist", ierr)
    close (NMLFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error('initialise :: init_namelist', "Failed to close "//trim(nmlfilename), ierr)

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
    FMT3, var2val(resuspension_coeff)
    FMT3, var2val(resuspension_threshold)
    FMT3, var2val(roughness_height)
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
    FMT3, var2val_char(xdimname)
    FMT3, var2val_char(ydimname)
    FMT3, var2val_char(zdimname)
    FMT3, var2val_char(uvarname)
    FMT3, var2val_char(vvarname)
    FMT3, var2val_char(zaxvarname)
    FMT3, var2val_char(elevvarname)
    FMT3, var2val_char(rhovarname)
    FMT3, var2val_char(tempvarname)
    FMT3, var2val_char(saltvarname)
    FMT3, var2val_char(viscvarname)
    FMT3, var2val_char(taubxvarname)
    FMT3, var2val_char(taubyvarname)
    FMT3, var2val(zax_style)
    FMT3, var2val(zax_direction)

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
    ! Initialise the fieldset
    ! TODO: Fix memory calculation
    !---------------------------------------------
    character(len=LEN_CHAR_L)              :: initPath
    character(len=LEN_CHAR_L)              :: filename
    character(len=LEN_CHAR_S), allocatable :: dimnames(:)
    real(rk)                               :: real_var
    ! integer(rk)                            :: field_count = 0, field_mem
    integer                                :: ndim      ! The (default) number of dimensions for the fieldset
    integer, allocatable                   :: dim_idx(:) ! The default dimensions for the fieldset

    FMT1, "======== Init fields ========"

    if (run_3d) then
      ndim = 3
      dim_idx = [1, 2, 3]
      dimnames = [character(LEN_CHAR_S) :: trim(xdimname), trim(ydimname), trim(zdimname)]
    else
      ndim = 2
      dim_idx = [1, 2]
      if (trim(zdimname) == "") then
        dimnames = [character(LEN_CHAR_S) :: trim(xdimname), trim(ydimname)]
      else
        dimnames = [character(LEN_CHAR_S) :: trim(xdimname), trim(ydimname), trim(zdimname)]
      end if
    end if

    ! dimnames = [character(LEN_CHAR_S) :: trim(xdimname), trim(ydimname), trim(zdimname)]

    ! field_mem = (nx * ny * nlevels * 2) * sizeof(real_var)
    ! FMT2, "Allocating fields of size (nx, ny, nz): (", nx, ", ", ny, ", ", nlevels, ")", field_mem, " bytes per field"
    if (has_subdomains) then
      fieldset = t_fieldset(path=GETMPATH, &
                            domain=domain, &
                            dimnames=dimnames, &
                            start=run_start_dt, dt=dt, &
                            file_prefix=trim(file_prefix), file_suffix=trim(file_suffix), &
                            pmap=PMAPFILE)
    else
      fieldset = t_fieldset(path=GETMPATH, &
                            domain=domain, &
                            dimnames=dimnames, &
                            start=run_start_dt, dt=dt, &
                            file_prefix=trim(file_prefix), file_suffix=trim(file_suffix))
    end if

    ! call fieldset%set_start_time(run_start_dt)
    ! call fieldset%set_simulation_timestep(dt)

    select case (has_subdomains)
    case (.true.)
      initPath = fieldset%get_directory(1)
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
    ! ! field_count = field_count + 2
    !---------------------------------------------
    ! Fields for v. velocity
    if (run_3d) then
      if (nc_var_exists(trim(filename), trim(wvarname))) then
        call fieldset%add_field("W", wvarname)
      else
        call throw_error("initialise :: init_fieldset", "Variable for W velocity does not exist: "//trim(wvarname))
      end if

      if (nc_var_exists(trim(filename), trim(zaxvarname))) then
        select case (zax_style)
        case (STATIC_DEPTH_VALUES)
          call fieldset%add_field("ZAX", zaxvarname)
        case (DEPTH_VALUES, LAYER_THICKNESS)
          call fieldset%add_field("ZAX", zaxvarname)
        end select
        call fieldset%set_zax("ZAX", zax_style, zax_direction)
      else
        call throw_error("initialise :: init_fieldset", "Variable for Z axis does not exist: "//trim(zaxvarname))
      end if

      ! ! field_count = field_count + 2

      if (nc_var_exists(trim(filename), trim(elevvarname))) then
        call fieldset%add_field("ELEV", elevvarname)
        ! TODO: set_elev for faster lookup
        ! ! field_count = field_count + 1
      end if
    end if
    !---------------------------------------------
    if (do_velocity) then
      !---------------------------------------------
      ! Field for density
      if (nc_var_exists(trim(filename), trim(rhovarname))) then
        call fieldset%add_field("RHO", rhovarname)
        ! ! field_count = field_count + 1
        density_method = RHO_VARIABLE
      else
        call throw_warning("initialise :: init_fieldset", "Could not find density ('"//trim(rhovarname)//"') in "//trim(filename))
        if (nc_var_exists(trim(filename), trim(tempvarname)) .and. &
            nc_var_exists(trim(filename), trim(saltvarname))) then
          call fieldset%add_field("TEMP", tempvarname)
          call fieldset%add_field("SALT", saltvarname)
          ! ! field_count = field_count + 2
          density_method = RHO_CALC
        else
          call throw_warning("initialise :: init_fieldset", "Could not find temperature or salinity ('" &
                             //trim(tempvarname)//"'/'"//trim(saltvarname)//"') in "//trim(filename)// &
                             ". Using default density.")
          density_method = RHO_DEFAULT
        end if
      end if
      !---------------------------------------------
      ! Field for viscosity
      if (nc_var_exists(trim(filename), trim(viscvarname))) then
        call fieldset%add_field("VISC", viscvarname)
        viscosity_method = VISC_VARIABLE
        ! ! field_count = field_count + 1
      else
       call throw_warning("initialise :: init_fieldset", "Could not find viscosity ('"//trim(viscvarname)//"') in "//trim(filename))
        if (nc_var_exists(trim(filename), trim(tempvarname)) .and. &
            nc_var_exists(trim(filename), trim(saltvarname))) then
          if (.not. fieldset%has_field("TEMP")) then
            call fieldset%add_field("TEMP", tempvarname)
            ! ! field_count = field_count + 1
          end if
          if (.not. fieldset%has_field("SALT")) then
            call fieldset%add_field("SALT", saltvarname)
            ! ! field_count = field_count + 1
          end if
          viscosity_method = VISC_CALC
        else
          call throw_warning("initialise :: init_fieldset", "Could not find temperature or salinity ('" &
                             //trim(tempvarname)//"'/'"//trim(saltvarname)//"') in "//trim(filename)// &
                             ". Using default viscosity.")
          viscosity_method = VISC_DEFAULT
        end if
      end if
      !---------------------------------------------
      ! Fields for bottom friction
      if ((nc_var_exists(trim(filename), trim(taubxvarname)) .and. (nc_var_exists(trim(filename), trim(taubyvarname))))) then
        call fieldset%add_field("TAUBX", taubxvarname)
        call fieldset%add_field("TAUBY", taubyvarname)
        ! ! field_count = field_count + 2
        has_bottom_stress = .true.
      else
        call throw_warning("initialise :: init_fieldset", "Could not find bottom stress ('"//trim(taubxvarname)//"'/'"//trim(taubyvarname)//"') in "//trim(filename))
      end if
    end if

    !---------------------------------------------
    if (do_biofouling) then
      call init_biofouling(fieldset)
    end if
    !---------------------------------------------

    FMT1, LINE
    ! FMT2, "Fields allocated: total", fieldset%num_fields * field_mem, " bytes"
    call fieldset%list_fields()

  end subroutine init_fieldset
  !===========================================
  subroutine init_particles

    call initialise_particles(fieldset)

  end subroutine init_particles
  !===========================================
  subroutine init_model
    !---------------------------------------------
    ! Call all the subroutines to initialise the model
    !---------------------------------------------

    call init_namelist                 ! init.f90
    call init_domain                   ! init.f90
    call init_time                     ! init.f90
    call init_fieldset                 ! init.f90
    call init_particles                ! init.f90
    call init_output                   ! output.f90

  end subroutine init_model

end module mod_initialise
