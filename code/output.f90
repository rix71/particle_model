#include "cppdefs.h"
#include "file.h"
#include "ompdefs.h"
#include "output.h"
module mod_output
  !----------------------------------------------------------------
  ! Module for writing the output
  !----------------------------------------------------------------
#ifdef USE_OMP
  use omp_lib
#endif
  use mod_precdefs
  use mod_errors
  use mod_particle
  use run_params, only: runid, nmlfilename
  use mod_params, only: run_3d
  use mod_particle_vars, only: particles, n_total_particles
  use time_vars, only: theDate, nTimes, dt
  use netcdf
  use nc_manager, only: nc_write, nc_initialise, &
                        nc_add_dimension, nc_add_variable, &
                        nc_add_attr, nc_check, &
                        FILLVALUE_BIG
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: outputstep, restartstep, init_output, &
            write_data, &
            write_restart, &
            write_all_particles, &
            write_active_particles, &
            outDir
#ifdef OUT_SUPPORT_ACTIVE
  public :: write_data_only_active
#endif
  !---------------------------------------------
  integer                   :: outputstep, restartstep
  character(len=LEN_CHAR_L) :: outDir
  integer                   :: nc_t_dimid, nc_p_dimid
  character(len=LEN_CHAR_L) :: nc_fileout_all
#ifdef OUT_SUPPORT_ACTIVE
  character(len=LEN_CHAR_L) :: nc_fileout_active
#endif
  logical                   :: write_all_particles, write_active_particles
  namelist /output_vars/ outDir, outputstep, restartstep, write_all_particles, write_active_particles
  !---------------------------------------------
  integer                     :: ierr
  !===================================================
contains
  !===========================================
  subroutine init_output
    logical :: dirExists

    FMT1, "======== Init output ========"

    open (NMLFILE, file=trim(nmlfilename), action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error("output :: init_output", "Failed to open "//trim(nmlfilename), ierr)
    read (NMLFILE, nml=output_vars, iostat=ierr)
    if (ierr .ne. 0) call throw_error("output :: init_output", "Failed to read "//trim(nmlfilename), ierr)
    close (NMLFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error("output :: init_output", "Failed to close "//trim(nmlfilename), ierr)

    inquire (file=trim(outDir), exist=dirExists)
    if (.not. dirExists) then
#ifndef NOSYSCALLS
      FMT2, "Making directory "//trim(outDir)
      call system('mkdir -p '//trim(outDir))
#else
      call throw_error("output :: init_output", "Out dir ("//trim(outDir)//") does not exist!")
#endif
    end if

    FMT2, var2val(write_all_particles)
    FMT2, var2val(write_active_particles)
    FMT2, "Writing output every ", outputstep, " timesteps, or ", (outputstep * dt) / 3600., "hours"
    FMT2, "Saving ", nTimes / outputstep, " timesteps"
    if (restartstep > 0) then
      FMT2, "Writing restart every ", restartstep, " timesteps, or ", (restartstep * dt) / 3600., "hours"
    else if (restartstep == 0) then
      FMT2, "Writing restart at end of simulation"
    end if

    if (write_all_particles) then
      nc_fileout_all = trim(outDir)//'/'//trim(runid)//'.out.nc'
      call init_nc_output(nc_fileout_all)
    end if

    if (write_active_particles) then
#ifdef OUT_SUPPORT_ACTIVE
      nc_fileout_active = trim(outDir)//'/'//trim(runid)//'.active.out.nc'
      call init_nc_output(nc_fileout_active)
#else
      call throw_warning("output :: init_output", "write_active_particles is true, but OUT_SUPPORT_ACTIVE is not defined")
#endif
    end if

    return
  end subroutine init_output
  !===========================================
  subroutine init_nc_output(file_name)

    character(len=LEN_CHAR_L), intent(in) :: file_name

    call nc_initialise(file_name)
    call nc_add_dimension(file_name, "time", nc_t_dimid)
    call nc_add_dimension(file_name, "particle", nc_p_dimid, n_total_particles)

    call nc_add_variable(file_name, "time", "float", 1, [nc_t_dimid])
    call nc_add_attr(file_name, "time", "units", "seconds since 1900-01-01 00:00:00")

    call nc_add_variable(file_name, "lon", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "lon", "units", "degrees east")

    call nc_add_variable(file_name, "lat", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "lat", "units", "degrees north")

    call nc_add_variable(file_name, "depth", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "depth", "units", "m")
    call nc_add_attr(file_name, "depth", "name", "depth")
#ifdef OUT_ID
    call nc_add_variable(file_name, "id", "float", 1, [nc_p_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "id", "name", "particle id")
#endif
#ifdef OUT_VELOCITY
    call nc_add_variable(file_name, "vx", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "vx", "units", "m/s")
    call nc_add_attr(file_name, "vx", "name", "eastward velocity")

    call nc_add_variable(file_name, "vy", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "vy", "units", "m/s")
    call nc_add_attr(file_name, "vy", "name", "northward velocity")

    call nc_add_variable(file_name, "vz", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "vz", "units", "m/s")
    call nc_add_attr(file_name, "vz", "name", "vertical velocity")
#endif
#ifdef OUT_SETTLING_VELOCITY
    call nc_add_variable(file_name, "vs", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "vs", "units", "m/s")
    call nc_add_attr(file_name, "vs", "name", "vertical velocity (settling) calculated from density difference")
#endif
#ifdef OUT_DENSITY
    call nc_add_variable(file_name, "rho", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "rho", "units", "kg/m3")
    call nc_add_attr(file_name, "rho", "name", "density")
#endif
#ifdef OUT_DENSITY_PLASTIC
    call nc_add_variable(file_name, "rho_plastic", "float", 1, [nc_p_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "rho_plastic", "units", "kg/m3")
    call nc_add_attr(file_name, "rho_plastic", "name", "initial density of plastic particle")
#endif
#ifdef OUT_DELTA_RHO
    call nc_add_variable(file_name, "delta_rho", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "delta_rho", "units", "kg/m3")
    call nc_add_attr(file_name, "delta_rho", "name", "density difference")
#endif
#ifdef OUT_RADIUS
    call nc_add_variable(file_name, "radius", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "radius", "units", "m")
    call nc_add_attr(file_name, "radius", "name", "radius")
#endif
#ifdef OUT_RADIUS_PLASTIC
    call nc_add_variable(file_name, "radius_plastic", "float", 1, [nc_p_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "radius_plastic", "units", "m")
    call nc_add_attr(file_name, "radius_plastic", "name", "initial radius of plastic particle")
#endif
#ifdef OUT_AGE
    call nc_add_variable(file_name, "age", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "age", "units", "s")
    call nc_add_attr(file_name, "age", "name", "age of the particle")
#endif
#ifdef OUT_TRAJECTORY
    call nc_add_variable(file_name, "trajectory", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "trajectory", "units", "m")
    call nc_add_attr(file_name, "trajectory", "name", "distance travelled")
#endif
#ifdef OUT_TIME_ON_BEACH
    call nc_add_variable(file_name, "time_on_beach", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "time_on_beach", "units", "s")
    call nc_add_attr(file_name, "time_on_beach", "name", "time on beach")
#endif
#ifdef OUT_KIN_VISCOSITY
    call nc_add_variable(file_name, "kin_viscosity", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "kin_viscosity", "units", "m2/s")
    call nc_add_attr(file_name, "kin_viscosity", "name", "kinematic viscosity of the surrounding seawater")
#endif
#ifdef OUT_FRICTION_VELOCITY
    call nc_add_variable(file_name, "u_star", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "u_star", "units", "m/s")
    call nc_add_attr(file_name, "u_star", "name", "friction velocity")
#endif
#ifdef OUT_H_BIOFILM
    call nc_add_variable(file_name, "h_biofilm", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "h_biofilm", "units", "m")
    call nc_add_attr(file_name, "h_biofilm", "name", "thickness of biofilm")
#endif
#ifdef OUT_GROWTH_BIOFILM
    call nc_add_variable(file_name, "growth_biofilm", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "growth_biofilm", "units", "m/s")
    call nc_add_attr(file_name, "growth_biofilm", "name", "growth rate of biofilm")
#endif
#ifdef OUT_STATE
    call nc_add_variable(file_name, "state", "int", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "state", "name", "state of the particle: 1-beached, 2-on boundary, 3-active, 4-bottom")
#endif

  end subroutine init_nc_output
  !===========================================
  subroutine write_data(nwrite)
    !---------------------------------------------
    ! Write the output
    ! Gathers all data in a matrix and writes it in one go
    !---------------------------------------------
    integer, intent(in)   :: nwrite
    integer               :: ipart, ncid, varid
    integer, save         :: nc_itime_out = 0
    real(rk), allocatable :: var(:, :)
    integer, allocatable  :: var_int(:, :)
    real(rk)              :: dateval(1)

    call theDate%print_short_date
    FMT2, "Saving all... ", nwrite, " particles"

    nc_itime_out = nc_itime_out + 1
    call nc_check(trim(nc_fileout_all), nf90_open(trim(nc_fileout_all), nf90_write, ncid), "write_data :: open")
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "time", varid), "write_data :: inq_varid")
    dateval(1) = theDate%date2num()
   call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, dateval, start=[nc_itime_out], count=[1]), "write_data :: put_var")

#ifdef OUT_VELOCITY
    allocate (var(nwrite, 6))
#else
    allocate (var(nwrite, 3))
#endif

    GATHER_OMP_START_RV
    do ipart = 1, nwrite
      var(ipart, 1) = particles(ipart)%lon0
      var(ipart, 2) = particles(ipart)%lat0
      var(ipart, 3) = particles(ipart)%depth0
#ifdef OUT_VELOCITY
      var(ipart, 4) = particles(ipart)%u0
      var(ipart, 5) = particles(ipart)%v0
      var(ipart, 6) = particles(ipart)%w0
#endif
    end do
    GATHER_OMP_END

    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "lon", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,1), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "lat", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,2), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "depth", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,3), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
#ifdef OUT_VELOCITY
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "vx", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,4), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "vy", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,5), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "vz", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,6), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
#endif

    ! Deallocate the matrix for single variables
    deallocate (var)

#ifdef OUT_ID
    if (.not. allocated(var)) allocate (var(nwrite, 1))
    GATHER_OMP_START_RV
    GATHER_LOOP_RV(id)
    GATHER_OMP_END
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "id", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,1), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
    ! Set to zero to avoid deallocation
    ! deallocate (var)
    var = ZERO
#endif

#ifdef OUT_SETTLING_VELOCITY
    if (.not. allocated(var)) allocate (var(nwrite, 1))
    GATHER_OMP_START_RV
    GATHER_LOOP_RV(vel_vertical)
    GATHER_OMP_END
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "vs", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,1), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
    ! Set to zero to avoid deallocation
    ! deallocate (var)
    var = ZERO
#endif
#ifdef OUT_DENSITY
    if (.not. allocated(var)) allocate (var(nwrite, 1))
    GATHER_OMP_START_RV
    GATHER_LOOP_RV(rho)
    GATHER_OMP_END
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "rho", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,1), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
    ! Set to zero to avoid deallocation
    ! deallocate (var)
    var = ZERO
#endif
#ifdef OUT_DENSITY_PLASTIC
    if (.not. allocated(var)) allocate (var(nwrite, 1))
    GATHER_OMP_START_RV
    GATHER_LOOP_RV(rho0)
    GATHER_OMP_END
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "rho_plastic", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,1), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
    ! Set to zero to avoid deallocation
    ! deallocate (var)
    var = ZERO
#endif
#ifdef OUT_DELTA_RHO
    if (.not. allocated(var)) allocate (var(nwrite, 1))
    GATHER_OMP_START_RV
    GATHER_LOOP_RV(delta_rho)
    GATHER_OMP_END
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "delta_rho", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,1), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
    ! Set to zero to avoid deallocation
    ! deallocate (var)
    var = ZERO
#endif
#ifdef OUT_RADIUS
    if (.not. allocated(var)) allocate (var(nwrite, 1))
    GATHER_OMP_START_RV
    GATHER_LOOP_RV(radius)
    GATHER_OMP_END
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "radius", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,1), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
    ! Set to zero to avoid deallocation
    ! deallocate (var)
    var = ZERO
#endif
#ifdef OUT_RADIUS_PLASTIC
    if (.not. allocated(var)) allocate (var(nwrite, 1))
    GATHER_OMP_START_RV
    GATHER_LOOP_RV(radius0)
    GATHER_OMP_END
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "radius_plastic", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,1), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
    ! Set to zero to avoid deallocation
    ! deallocate (var)
    var = ZERO
#endif
#ifdef OUT_AGE
    if (.not. allocated(var)) allocate (var(nwrite, 1))
    GATHER_OMP_START_RV
    GATHER_LOOP_RV(age)
    GATHER_OMP_END
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "age", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,1), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
    ! Set to zero to avoid deallocation
    ! deallocate (var)
    var = ZERO
#endif
#ifdef OUT_TRAJECTORY
    if (.not. allocated(var)) allocate (var(nwrite, 1))
    GATHER_OMP_START_RV
    GATHER_LOOP_RV(traj_len)
    GATHER_OMP_END
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "trajectory", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,1), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
    ! Set to zero to avoid deallocation
    ! deallocate (var)
    var = ZERO
#endif
#ifdef OUT_TIME_ON_BEACH
    if (.not. allocated(var)) allocate (var(nwrite, 1))
    GATHER_OMP_START_RV
    GATHER_LOOP_RV(time_on_beach)
    GATHER_OMP_END
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "time_on_beach", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,1), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
    ! Set to zero to avoid deallocation
    ! deallocate (var)
    var = ZERO
#endif
#ifdef OUT_H_BIOFILM
    if (.not. allocated(var)) allocate (var(nwrite, 1))
    GATHER_OMP_START_RV
    GATHER_LOOP_RV(h_biofilm)
    GATHER_OMP_END
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "h_biofilm", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,1), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
    ! Set to zero to avoid deallocation
    ! deallocate (var)
    var = ZERO
#endif
#ifdef OUT_GROWTH_BIOFILM
    if (.not. allocated(var)) allocate (var(nwrite, 1))
    GATHER_OMP_START_RV
    GATHER_LOOP_RV(growth_biofilm)
    GATHER_OMP_END
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "growth_biofilm", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,1), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
    ! Set to zero to avoid deallocation
    ! deallocate (var)
    var = ZERO
#endif
#ifdef OUT_STATE
    if (.not. allocated(var_int)) allocate (var_int(nwrite, 1))
    GATHER_OMP_START_IV
    GATHER_LOOP_IV(state)
    GATHER_OMP_END
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "state", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var_int(:,1), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
    ! Set to zero to avoid deallocation
    ! deallocate (var_int)
    var_int = 0
#endif
#ifdef OUT_KIN_VISCOSITY
    if (.not. allocated(var)) allocate (var(nwrite, 1))
    GATHER_OMP_START_RV
    GATHER_LOOP_RV(kin_viscosity)
    GATHER_OMP_END
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "kin_viscosity", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,1), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
    ! Set to zero to avoid deallocation
    ! deallocate (var)
    var = ZERO
#endif
#ifdef OUT_FRICTION_VELOCITY
    if (.not. allocated(var)) allocate (var(nwrite, 1))
    GATHER_OMP_START_RV
    GATHER_LOOP_RV(u_star)
    GATHER_OMP_END
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "u_star", varid), "write_data :: inq_varid")
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, var(:,1), start=[1, nc_itime_out], count=[nwrite, 1]), "write_data :: put_var")
    ! Set to zero to avoid deallocation
    ! deallocate (var)
    var = ZERO
#endif

    call nc_check(trim(nc_fileout_all), nf90_close(ncid), "write_data :: close")

    deallocate (var)

    return
  end subroutine write_data
!===========================================
#ifdef OUT_SUPPORT_ACTIVE
  subroutine write_data_only_active(nwrite)

    integer, intent(in) :: nwrite
    integer             :: ipart, ncid, varid
    integer, save       :: nc_itime_out = 0
    real(rk)            :: var1d(1), var2d(1, 1), dateval(1)
    integer             :: var2d_int(1, 1)

    call theDate%print_short_date
    FMT2, "Saving active... ", nwrite, " particles"

    nc_itime_out = nc_itime_out + 1
    call nc_check(trim(nc_fileout_active), nf90_open(trim(nc_fileout_active), nf90_write, ncid), "write_data_active :: open")
    call nc_check(trim(nc_fileout_active), nf90_inq_varid(ncid, "time", varid), "write_data_active :: inq varid")
    dateval = theDate%date2num()
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, dateval, start=[nc_itime_out], count=[1]), &
                  "write_data_active :: put var")

    do ipart = 1, nwrite

      if (particles(ipart)%is_active) then

        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "lon", varid), "write_data :: inq varid")
        var2d = particles(ipart)%lon0
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data :: put var")

        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "lat", varid), "write_data :: inq varid")
        var2d = particles(ipart)%lat0
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data :: put var")

        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "depth", varid), "write_data :: inq varid")
        var2d = particles(ipart)%depth0
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data :: put var")
#ifdef OUT_ID
        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "id", varid), "write_data :: inq varid")
        var1d = particles(ipart)%id
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var1d, start=[ipart], count=[1]), &
                      "write_data :: put var")
#endif
#ifdef OUT_VELOCITY
        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "vx", varid), "write_data :: inq varid")
        var2d = particles(ipart)%u0
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data :: put var")

        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "vy", varid), "write_data :: inq varid")
        var2d = particles(ipart)%v0
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data :: put var")

        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "vz", varid), "write_data :: inq varid")
        var2d = particles(ipart)%w0
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data :: put var")
#endif
#ifdef OUT_SETTLING_VELOCITY
        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "vs", varid), "write_data :: inq varid")
        var2d = particles(ipart)%vel_vertical
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data :: put var")
#endif
#ifdef OUT_DENSITY
        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "rho", varid), "write_data :: inq varid")
        var2d = particles(ipart)%rho
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data :: put var")
#endif
#ifdef OUT_DENSITY_PLASTIC
        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "rho_plastic", varid), "write_data :: inq varid")
        var1d = particles(ipart)%rho0
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart], count=[1]), &
                      "write_data :: put var")
#endif
#ifdef OUT_DELTA_RHO
        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "delta_rho", varid), "write_data :: inq varid")
        var2d = particles(ipart)%delta_rho
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data :: put var")
#endif
#ifdef OUT_RADIUS
        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "radius", varid), "write_data :: inq varid")
        var2d = particles(ipart)%radius
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data :: put var")
#endif
#ifdef OUT_RADIUS_PLASTIC
        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "radius_plastic", varid), "write_data :: inq varid")
        var1d = particles(ipart)%radius0
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var1d, start=[ipart], count=[1]), &
                      "write_data :: put var")
#endif
#ifdef OUT_AGE
        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "age", varid), "write_data :: inq varid")
        var2d = particles(ipart)%age
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data :: put var")
#endif
#ifdef OUT_TRAJECTORY
        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "trajectory", varid), "write_data :: inq varid")
        var2d = particles(ipart)%traj_len
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data :: put var")
#endif
#ifdef OUT_TIME_ON_BEACH
        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "time_on_beach", varid), "write_data :: inq varid")
        var1d = particles(ipart)%time_on_beach
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart], count=[1]), &
                      "write_data :: put var")
#endif
#ifdef OUT_KIN_VISCOSITY
        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "kin_viscosity", varid), "write_data :: inq varid")
        var2d = particles(ipart)%kin_visc
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data :: put var")
#endif
#ifdef OUT_FRICTION_VELOCITY
        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "u_star", varid), "write_data :: inq varid")
        var2d = particles(ipart)%u_star
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data :: put var")
#endif
#ifdef OUT_H_BIOFILM
        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "h_biofilm", varid), "write_data :: inq varid")
        var2d = particles(ipart)%h_biofilm
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data :: put var")
#endif
#ifdef OUT_GROWTH_BIOFILM
        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "growth_biofilm", varid), "write_data :: inq varid")
        var2d = particles(ipart)%growth_biofilm
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data :: put var")
#endif
#ifdef OUT_STATE
        call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "state", varid), "write_data :: inq varid")
        var2d_int = particles(ipart)%state
        call nc_check(trim(nc_fileout_all), &
                      nf90_put_var(ncid, varid, var2d_int, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data :: put var")
#endif
      end if
    end do

    call nc_check(trim(nc_fileout_active), nf90_close(ncid), "write_data_active :: close")

  end subroutine write_data_only_active
#endif
  !===========================================
  subroutine write_restart(nwrite)
    ! ---------------------------------------
    ! Write restart file
    ! TODO: this should be a method of the particle (array) class
    ! ---------------------------------------
    integer, intent(in) :: nwrite
    character(len=LEN_CHAR_L) :: restart_file
    character(len=14) :: time_str
    integer :: ipart

    call theDate%print_short_date
    FMT2, "Saving restart... ", nwrite, " particles"

    write (time_str, '(i0.14)') theDate%shortDate(include_time=.true.)
    restart_file = trim(outDir)//"/"//trim(runid)//"."//trim(time_str)//".restart.dat"

    open (RESTARTFILE, file=trim(restart_file), action='write', status='new', iostat=ierr)
    write (RESTARTFILE, *) nwrite

    do ipart = 1, nwrite
      write (RESTARTFILE, *) particles(ipart)%lon0, particles(ipart)%lat0, particles(ipart)%depth0, &
        particles(ipart)%i0, particles(ipart)%j0, particles(ipart)%k0, &
        particles(ipart)%ir0, particles(ipart)%jr0, particles(ipart)%kr0, &
        particles(ipart)%id, particles(ipart)%beaching_time, &
        particles(ipart)%rho, particles(ipart)%rho0, &
        particles(ipart)%radius, particles(ipart)%radius0, &
        particles(ipart)%h_biofilm, &
        particles(ipart)%age, particles(ipart)%max_age, particles(ipart)%kill_beached, particles(ipart)%kill_boundary, &
        particles(ipart)%u0, particles(ipart)%v0, particles(ipart)%w0, particles(ipart)%vel_vertical, &
        particles(ipart)%traj_len, particles(ipart)%time_on_beach, particles(ipart)%is_active, particles(ipart)%state
    end do
    close (RESTARTFILE)

  end subroutine write_restart
end module mod_output
