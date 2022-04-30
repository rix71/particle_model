#include "cppdefs.h"
module mod_loop
  !----------------------------------------------------------------
  ! Main loop
  !----------------------------------------------------------------
#ifdef USE_OMP
  use omp_lib
#endif
  use mod_precdefs
  use mod_errors
  use mod_params, only: do_velocity, do_diffusion, run_3d, advection_method
  use mod_advection
  use mod_physics
  use field_vars, only: fieldset, has_density, has_viscosity
  use mod_particle, only: t_particle
  use mod_domain_vars, only: domain
  use mod_particle_vars, only: particles, init_coords, inputstep, &
                               max_age, runparts, kill_beached, kill_boundary
  use time_vars, only: theDate, run_start_dt, run_end_dt, dt
  use mod_output, only: outputstep, snap_interval, write_data, &
                        write_beached, write_boundary, write_data_only_active, &
                        write_all_particles, write_active_particles, write_data_snapshot, write_snapshot
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: loop
  !---------------------------------------------
  !===================================================
contains
  !===========================================
  subroutine loop
    integer           :: itime = 0
    integer           :: ipart, i_release = 1
    real(rk)          :: time
    character(len=8)  :: d
    character(len=10) :: t
#ifndef SAY_LESS
    integer           :: active_particles, inactive_particles
#endif

    dbghead(loop)

    FMT1, "======== Starting time loop ========"
    call date_and_time(date=d, time=t)
    FMT2, t(1:2), ":", t(3:4), ":", t(5:10)
#ifdef USE_OMP
    FMT2, "Using OpenMP with ", omp_get_num_procs(), " processes"
#endif

    ! Read appropriate fields:
    call fieldset%read_first_timesteps(run_start_dt)

    ! Start integration loop
    do while (theDate < run_end_dt)

      !   - update fields
      call fieldset%update(theDate)
      time = fieldset%get_time(theDate)
#ifdef SMAGORINSKY_FULL_FIELD
      if (do_diffusion) call update_Ah_Smagorinsky_full_field(fieldset, time, fieldset%nx, fieldset%ny, fieldset%nz)
#endif

      !   - release particles
      if (mod(itime, inputstep) .eq. 0) then
        FMT2, "Releasing ", init_coords(i_release)%n_particles, " new particles at itime = ", itime
        do ipart = 1, init_coords(i_release)%n_particles
          particles(ipart + runparts) = t_particle(lon=init_coords(i_release)%x(ipart), &
                                                   lat=init_coords(i_release)%y(ipart), &
                                                   depth=init_coords(i_release)%z(ipart), &
                                                   id=init_coords(i_release)%id(ipart), &
                                                   beaching_time=init_coords(i_release)%beaching_time(ipart), &
                                                   rho=init_coords(i_release)%rho(ipart), &
                                                   radius=init_coords(i_release)%radius(ipart), &
                                                   max_age=max_age, &
                                                   kill_beached=kill_beached, &
                                                   kill_boundary=kill_boundary, &
                                                   fieldset=fieldset, &
                                                   time=time)
        end do
        runparts = runparts + init_coords(i_release)%n_particles
        FMT2, runparts, "particles"
        i_release = init_coords(i_release)%next_idx
      end if

#ifndef SAY_LESS
      if (mod(itime, PROGRESSINFO) .eq. 0) then
        call theDate%print_short_date
        call date_and_time(date=d, time=t)
        FMT2, t(1:2), ":", t(3:4), ":", t(5:10), " itime = ", itime
        if (itime .ne. 0) then
          active_particles = 0
          inactive_particles = 0
#ifdef USE_OMP
          !$omp parallel do reduction(+: active_particles, inactive_particles)
#endif
          do ipart = 1, runparts
            if (particles(ipart)%is_active) then
              active_particles = active_particles + 1
            else
              inactive_particles = inactive_particles + 1
            end if
          end do
#ifdef USE_OMP
          !$omp end parallel do
#endif
          FMT2, active_particles, " active particles"
          FMT2, inactive_particles, " inactive particles"
        end if
      end if
#endif

      !---------------------------------------------
      ! Start particle loop
#ifdef USE_OMP
#ifdef SMAGORINSKY_FULL_FIELD
      !$omp parallel do shared(particles, fieldset, domain, Ah_field)
#else
      !$omp parallel do shared(particles, fieldset, domain)
#endif
#endif
      do ipart = 1, runparts
        DBG, "----------------------------------------"
        DBG, "    Particle nr: ", ipart, "            "
        DBG, "----------------------------------------"
        !---------------------------------------------
        ! Skip inactive particles
        if (.not. particles(ipart)%is_active) cycle
        DBG, "Still active :)"
        !---------------------------------------------
        ! Advect only if the particle is alive (is_active=.true.) and active (state=0)
        if (particles(ipart)%state == ACTIVE) then
          ! - do advection
          call advect(particles(ipart), fieldset, time, advection_method, run_3d)
          ! - do Kooi vertical velocity
          if (do_velocity .and. run_3d) call vertical_velocity(particles(ipart), fieldset, time, has_density, has_viscosity)
          ! - do diffusion
          if (do_diffusion) call diffuse(particles(ipart), fieldset, time, run_3d)
        end if

        !---------------------------------------------
        ! Update particles
#ifdef DEBUG
        DBG, "PARTICLE BEFORE UPDATING: ", ipart
        call particles(ipart)%print_info()
#endif

        call particles(ipart)%update(fieldset, time)

#ifdef DEBUG
        DBG, "PARTICLE AFTER UPDATING: ", ipart
        call particles(ipart)%print_info()
#endif

#ifndef USE_OMP
        !---------------------------------------------
        ! Write snapshot
        ! Cannot use this with openMP (parallel io?)
        if ((mod(particles(ipart)%age, snap_interval) == 0) .and. (write_snapshot)) then
          call write_data_snapshot(particles(ipart), ipart)
        end if
#endif
      end do
#ifdef USE_OMP
      !$omp end parallel do
#endif

      !---------------------------------------------
      ! Write output
      if ((mod(itime, outputstep) .eq. 0) .and. (runparts .gt. 0)) then
        if (write_all_particles) call write_data(runparts)
        if (write_active_particles) call write_data_only_active(runparts)
      end if
      !---------------------------------------------
      ! Update time
      call theDate%update(dt)
      itime = itime + 1
    end do

    call date_and_time(date=d, time=t)
    FMT2, LINE
    FMT2, t(1:2), ":", t(3:4), ":", t(5:10)
    FMT2, "Finished all time steps (", itime, ")"
    FMT2, LINE

    dbgtail(loop)
    return
  end subroutine loop
end module mod_loop
