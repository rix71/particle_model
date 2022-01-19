#include "cppdefs.h"
module loop_particle
  !----------------------------------------------------------------
  ! Main loop
  !----------------------------------------------------------------
  use precdefs
  use errors
  use interp
  use params, only: do_diffusion, do_velocity
  use physics, only: diffuse, velocity, vertical_velocity, seawater_density_from_temp_and_salt
  use loop_vars, only: dateThis, dateNext, dateThisNcTimestep, dateNextNcTimestep, thisPath, nextPath, &
                       xnow, xnew, ynow, ynew, znow, znew, cartx, cartx_new, &
                       carty, carty_new, ig, jg, kg, igr, jgr, kgr, &
                       nc_itime, nc_itime_next, ncNTimes, pvelu, pvelunew, &
                       pvelv, pvelvnew, pvelw, pvelwnew
  use particle_type, only: particle
  use particle_vars, only: nInitParticles, runparts, inputstep, initCoords, particles, max_age, kill_beached, kill_boundary
  use time_vars, only: theDate, run_end_dt, dt, nc_timestep, nTimes
  use domain, only: lonlat2xy, xy2lonlat
  use domain_vars, only: x0, y0, dx, dy, dz, dx_m, dy_m, nx, ny, lons, lats, depdata, seamask
  use field_vars, only: nlevels, has_subdomains, has_viscosity, has_density, run_3d, &
                        uspeed, vspeed, wspeed, &
                        uspeednew, vspeednew, wspeednew, &
                        udata, vdata, wdata, &
                        udatanew, vdatanew, wdatanew, &
                        udata_interp, vdata_interp, wdata_interp, &
                        udatanew_interp, vdatanew_interp, wdatanew_interp, &
                        zaxdata, zaxdatanew, zaxdata_interp, zaxdatanew_interp, &
                        density, densitynew, temp, tempnew, salt, saltnew, visc, viscnew, &
                        rho_sea, viscosity, file_prefix, file_suffix
  use fields, only: find_folder, find_file, get_indices2d, &
                    get_indices_vert, read_fields_full_domain
  use nc_manager, only: nc_get_dim
  use modtime
  use output, only: outputstep, write_data, &
                    write_beached, write_boundary, write_data_only_active, &
                    write_all_particles, write_active_particles, write_data_snapshot
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: loop
  !---------------------------------------------
  integer           :: active_particles, inactive_particles
  character(len=8)  :: d
  character(len=10) :: t
  !===================================================
contains
  !===========================================
  subroutine loop
    logical :: read_first = .false., interp_first = .false. ! This should avoid reading the fields twice
    integer :: ipart, itime = 0, itime_interp = 0, n_adjust_dt = 0
    integer :: ig_prev, jg_prev
    real(rk) :: dt_orig, dt_interp_1, dt_interp_2

    dbghead(loop)

    call print_loop_start
    !---------------------------------------------
    ! Initialise netCDF first dates
    select case (has_subdomains)
    case (.true.)
      call get_date_and_path_dirlist
    case (.false.)
      call get_date_and_path_filelist
    end select
    !---------------------------------------------
    ! Save a copy of the original timestep
    dt_orig = dt
    !---------------------------------------------
    ! Start time loop
    !   The condition must be < and not <= because
    !   the date is updated in the interpolation loop
    !   and may result in an infinite loop in some cases.
    do while (theDate < run_end_dt)
      !---------------------------------------------
      ! Get a new path when it is time to move to the next folder.
      ! While you're at it, also update dateThis and dateNext.
      ! TODO: This should not be the only way to find the data,
      ! e.g. the data may be stored in one file or folder, or different names...
      if (theDate >= dateNext) then
        select case (has_subdomains)
        case (.true.)
          FMT2, "Moving to next folder..."
          call get_date_and_path_dirlist
        case (.false.)
          FMT2, "Moving to next file..."
          call get_date_and_path_filelist
        end select
      end if
      !---------------------------------------------
      ! Get netCDF time index
      debug(thisPath)
      debug(nextPath)
      select case (has_subdomains)
      case (.true.)
        if (dateThis == init_datetime_from_netcdf(trim(thisPath)//PROC0, 1)) then
          ! If RefTime is the same as the first timestep, then 1 must be added
          nc_itime = int(dateDiff(dateThis, theDate) / nc_timestep) + 1
        else
          nc_itime = int(dateDiff(dateThis, theDate) / nc_timestep)
        end if
      case (.false.)
        if (dateThis == init_datetime_from_netcdf(trim(thisPath), 1)) then
          ! If RefTime is the same as the first timestep, then 1 must be added
          nc_itime = int(dateDiff(dateThis, theDate) / nc_timestep) + 1
        else
          nc_itime = int(dateDiff(dateThis, theDate) / nc_timestep)
        end if
      end select

      nc_itime_next = nc_itime + 1
      debug(nc_itime)
      debug(nc_itime_next)
      debug(ncNTimes)
#ifdef DEBUG
      DBG, "dateThis:"; call dateThis%print_short_date
      DBG, "dateNext:"; call dateNext%print_short_date
#endif

      if (nc_itime_next .gt. ncNTimes) then
        DBG, "Updating path"
        !---------------------------------------------
        nc_itime_next = nc_itime_next - ncNTimes
        select case (has_subdomains)
        case (.true.)
          call find_folder(dateNext, nextPath)
        case (.false.)
          call find_file(dateNext, nextPath)
        end select
        debug(thisPath)
        debug(nextPath)
        debug(nc_itime)
        debug(nc_itime_next)
      end if

      !---------------------------------------------
      ! Date of the next timestep
      select case (has_subdomains)
      case (.true.)
        dateThisNcTimestep = init_datetime_from_netcdf(trim(thisPath)//PROC0, &
                                                       nc_itime)
        dateNextNcTimestep = init_datetime_from_netcdf(trim(nextPath)//PROC0, &
                                                       nc_itime_next)
      case (.false.)
        dateThisNcTimestep = init_datetime_from_netcdf(trim(thisPath), nc_itime)
        dateNextNcTimestep = init_datetime_from_netcdf(trim(nextPath), nc_itime_next)
      end select
#ifdef DEBUG
      DBG, "theDate:"; call theDate%print_short_date
      DBG, "dateThisNcTimestep:"; call dateThisNcTimestep%print_short_date
      DBG, "dateNextNcTimestep:"; call dateNextNcTimestep%print_short_date
#endif

      !---------------------------------------------
      ! Read in currents
      ! Does this even make sense???
      if (.not. read_first) then
        call read_fields_full_domain(nc_itime, thisPath, read_first)
        read_first = .true.
        call read_fields_full_domain(nc_itime_next, nextPath, read_first)
      else
        call read_fields_full_domain(nc_itime_next, nextPath, read_first)
      end if

      itime_interp = 0
      !---------------------------------------------
      ! Start interpolated time loop
      do while (theDate < dateNextNcTimestep)
        !---------------------------------------------
        if (itime >= nTimes) exit
        !---------------------------------------------
        ! Release particles
        ! TODO: Do this only in case of continuous release (particle_method ?)
        if (mod(itime, inputstep) .eq. 0) then
          FMT2, "Releasing new particles at itime = ", itime
          do ipart = 1, nInitParticles
            particles(ipart + runparts) = particle(xPos=initCoords(ipart, 1), &
                                                   yPos=initCoords(ipart, 2), &
                                                   originNum=initCoords(ipart, 3), &
                                                   beachingtime=initCoords(ipart, 4), &
                                                   rho=initCoords(ipart, 5), &
                                                   radius=initCoords(ipart, 6), &
                                                   kill_bch=kill_beached, &
                                                   kill_bdy=kill_boundary)
          end do
          runparts = runparts + nInitParticles
          FMT2, runparts, "particles"
        end if
        !---------------------------------------------
#ifndef SAYLESS
        if (mod(itime, PROGRESSINFO) .eq. 0) then
          call theDate%print_short_date
          call date_and_time(date=d, time=t)
          FMT2, t(1:2), ":", t(3:4), ":", t(5:10), " itime = ", itime
          if (itime .ne. 0) then
            active_particles = 0
            inactive_particles = 0
            do ipart = 1, runparts
              if (particles(ipart)%isActive) then
                active_particles = active_particles + 1
              else
                inactive_particles = inactive_particles + 1
              end if
            end do
            ! FMT2, "There are ", active_particles, " active particles"
            FMT2, active_particles, " active particles"
            FMT2, inactive_particles, " inactive particles"
          end if
        end if
#endif
        !---------------------------------------------
        ! Interpolate between timesteps
        if (dt .ne. dt_orig) then
          !---------------------------------------------
          ! If the next date is greater than the next nc date,
          ! then interpolate to exactly the next nc date and
          ! adjust the timestep
          ! Must go back to original timestep:
          !     dt_orig = dt_a + dt_b
          !     dt_1 = dt_a
          !     dt_2 = dt_b
          !     dt_3 = dt_orig
          ! This only happens if nc_timestep/dt_orig is not integer
          select case (n_adjust_dt)
          case (1)
            dt = dt_orig - dt
            n_adjust_dt = 2
          case (2)
            dt = dt_orig
            n_adjust_dt = 0
          end select
        end if
        dt_interp_1 = dateDiff(dateThisNcTimestep, theDate)
        dt_interp_2 = dt_interp_1 + dt
        if (dt_interp_2 > nc_timestep) then
          dt = nc_timestep - dt_interp_1
          n_adjust_dt = 1
        end if
        if (.not. interp_first) then
          call timeinterp(udata, udatanew, udata_interp, dt_interp_1, nc_timestep, nx, ny, nlevels, 1)
          call timeinterp(vdata, vdatanew, vdata_interp, dt_interp_1, nc_timestep, nx, ny, nlevels, 1)
          if (run_3d) then
            call timeinterp(wdata, wdatanew, wdata_interp, dt_interp_1, nc_timestep, nx, ny, nlevels, 1)
            call timeinterp(zaxdata, zaxdatanew, zaxdata_interp, dt_interp_1, nc_timestep, nx, ny, nlevels, 1)
          end if
          interp_first = .true.
        end if
        call timeinterp(udata, udatanew, udatanew_interp, min(dt_interp_2, nc_timestep), nc_timestep, nx, ny, nlevels, 1)
        call timeinterp(vdata, vdatanew, vdatanew_interp, min(dt_interp_2, nc_timestep), nc_timestep, nx, ny, nlevels, 1)
        if (run_3d) then
          call timeinterp(wdata, wdatanew, wdatanew_interp, min(dt_interp_2, nc_timestep), nc_timestep, nx, ny, nlevels, 1)
          call timeinterp(zaxdata, zaxdatanew, zaxdatanew_interp, min(dt_interp_2, nc_timestep), nc_timestep, nx, ny, nlevels, 1)
        end if
        !---------------------------------------------
        ! Start particle loop
        do ipart = 1, runparts
          DBG, "----------------------------------------"
          DBG, "    Particle nr: ", ipart, "            "
          DBG, "----------------------------------------"
          !---------------------------------------------
          ! Skip inactive particles
          if (.not. particles(ipart)%isActive) cycle
          DBG, "Still active :)"
          !---------------------------------------------
          ! Particle coordinates
          xnow = particles(ipart)%xPos; ynow = particles(ipart)%yPos
          znow = particles(ipart)%zPos
          !---------------------------------------------
          ! Advect only if the particle is alive (isActive=.true.) and active (state=0)
          if (particles(ipart)%state .eq. 0) then
            !---------------------------------------------
            ! Cartesian coordinates
            call lonlat2xy(xnow, ynow, x0, y0, cartx, carty)
            !---------------------------------------------
            ! Global i and j indices
            call get_indices2d(xnow, ynow, x0, y0, dx, dy, ig, jg)
            ig_prev = ig; jg_prev = jg
            !---------------------------------------------
            ! Get current speed at particle location
            select case (run_3d)
            case (.true.)
              call get_current_speed(xloc=xnow, yloc=ynow, zloc=znow, &
                                     uarr=udata_interp, varr=vdata_interp, warr=wdata_interp, zaxarr=zaxdata_interp, &
                                     uspeedout=uspeed, vspeedout=vspeed, wspeedout=wspeed)
              DBG, "Speeds 3D:"
              debug(uspeed); debug(vspeed); debug(wspeed)
            case (.false.)
              call get_current_speed(xloc=xnow, yloc=ynow, &
                                     uarr=udata_interp, varr=vdata_interp, &
                                     uspeedout=uspeed, vspeedout=vspeed)
              DBG, "Speeds 2D:"
              debug(uspeed); debug(vspeed)
            end select
            !---------------------------------------------
            ! Particle velocity
            select case (do_velocity)
            case (.true.)
              !---------------------------------------------
              ! Get seawater density
              select case (has_density)
              case (DEFAULT_DENSITY)
                ! Use some default value
                rho_sea = 1000.0d0
              case (DENSITY)
                rho_sea = density(ig, jg, kg)
              case (TEMP_SALT)
                rho_sea = seawater_density_from_temp_and_salt(temp(ig, jg, kg), salt(ig, jg, kg), znow)
              end select
              debug(rho_sea)
              !---------------------------------------------
              ! Get seawater viscosity
              select case (has_viscosity)
              case (.false.)
                ! Use some default value
              case (.true.)
                viscosity = visc(ig, jg, kg)
              end select
              debug(viscosity)
              !---------------------------------------------
         pvelu = velocity(xnow, ynow, particles(ipart)%u, particles(ipart)%rho, particles(ipart)%radius, uspeed, rho_sea, viscosity)
         pvelv = velocity(xnow, ynow, particles(ipart)%v, particles(ipart)%rho, particles(ipart)%radius, vspeed, rho_sea, viscosity)
              if (run_3d) pvelw = vertical_velocity(particles(ipart)%rho, particles(ipart)%radius, rho_sea)
            case (.false.)
              pvelu = uspeed
              pvelv = vspeed
              pvelw = wspeed
            end select
            DBG, "Particle velocity:"
            debug(pvelu); debug(pvelv); debug(pvelw)
            !---------------------------------------------
            ! Advect particle (first)
            cartx_new = cartx + pvelu * dt
            carty_new = carty + pvelv * dt
            if (run_3d) then
              znew = znow + pvelw * dt 
            end if
            !---------------------------------------------
            ! New coordinates and indices
            call xy2lonlat(cartx_new, carty_new, x0, y0, xnew, ynew)
            debug(xnow); debug(ynow); debug(znow)
            debug(xnew); debug(ynew); debug(znew)
            call get_indices2d(xnew, ynew, x0, y0, dx, dy, ig, jg)
            debug(ig); debug(jg)
            if (ig .gt. nx) xnew = lons(nx) ! TODO: These checks should go somewhere else
            if (jg .gt. ny) ynew = lats(ny) !       Is this even legal?
            !---------------------------------------------
            ! Particle should not be on land
            if (depdata(ig, jg) .lt. 0.0d0) then
              !---------------------------------------------
              ! Trying a different approach
              DBG, "Bouncing"
              debug(seamask(ig_prev, jg_prev))
              debug(seamask(ig, jg))
              if (ig .ne. ig_prev) then
                ! Particle has moved to the left or right box
                cartx_new = cartx - pvelu * dt
              end if
              if (jg .ne. jg_prev) then
                ! Particle has moved up or down
                carty_new = carty - pvelv * dt
              end if
              call xy2lonlat(cartx_new, carty_new, x0, y0, xnew, ynew)
              debug(xnow); debug(ynow); debug(znow)
              debug(xnew); debug(ynew); debug(znew)
              call get_indices2d(xnew, ynew, x0, y0, dx, dy, ig, jg)
              ig_prev = ig; jg_prev = jg
            end if
            !---------------------------------------------
            call particles(ipart)%check_depth(zaxdatanew_interp, ig, jg, znew, .false.)
            !---------------------------------------------
            ! Get current speed at particle location
            select case (run_3d)
            case (.true.)
              call get_current_speed(xloc=xnew, yloc=ynew, zloc=znew, &
                                     uarr=udatanew_interp, varr=vdatanew, warr=wdatanew_interp, zaxarr=zaxdatanew_interp, &
                                     uspeedout=uspeednew, vspeedout=vspeednew, wspeedout=wspeednew)
              DBG, "Speeds 3D new:"
              debug(uspeednew); debug(vspeednew); debug(wspeednew)
            case (.false.)
              call get_current_speed(xloc=xnew, yloc=ynew, &
                                     uarr=udatanew_interp, varr=vdatanew, &
                                     uspeedout=uspeednew, vspeedout=vspeednew)
              DBG, "Speeds 2D new:"
              debug(uspeednew); debug(vspeednew)

            end select
            !---------------------------------------------
            ! Particle velocity
            select case (do_velocity)
            case (.true.)
              !---------------------------------------------
              ! Get seawater density
              select case (has_density)
              case (DEFAULT_DENSITY)
                ! Use some default value
                rho_sea = 1000.0d0
              case (DENSITY)
                rho_sea = densitynew(ig, jg, kg)
              case (TEMP_SALT)
                rho_sea = seawater_density_from_temp_and_salt(tempnew(ig, jg, kg), saltnew(ig, jg, kg), znew)
              end select
              !---------------------------------------------
              ! Get seawater viscosity
              select case (has_viscosity)
              case (.false.)
                ! Use some default value
              case (.true.)
                viscosity = viscnew(ig, jg, kg)
              end select
              debug(viscosity)
              !---------------------------------------------
              pvelunew = velocity(xnew, ynew, pvelu, particles(ipart)%rho, particles(ipart)%radius, uspeednew, rho_sea, viscosity)
              pvelvnew = velocity(xnew, ynew, pvelv, particles(ipart)%rho, particles(ipart)%radius, vspeednew, rho_sea, viscosity)
              if (run_3d) pvelwnew = vertical_velocity(particles(ipart)%rho, particles(ipart)%radius, rho_sea)
            case (.false.)
              pvelunew = uspeednew
              pvelvnew = vspeednew
              pvelwnew = wspeednew
            end select
            DBG, "Particle velocity new:"
            debug(pvelunew); debug(pvelvnew); debug(pvelwnew)
            !---------------------------------------------
            ! Advect particle (second)
            cartx_new = cartx + 0.5 * (pvelu + pvelunew) * dt
            carty_new = carty + 0.5 * (pvelv + pvelvnew) * dt
            if (run_3d) then
              znew = znow + 0.5 * (pvelw + pvelwnew) * dt ! TODO: Don't let the particle jump out of water
              if (znew .gt. 0.0d0) then
                DBG, "Setting znew to 0 again"; debug(znew)
                znew = 0.0d0
              end if
            end if
            !---------------------------------------------
            ! Particle should not be on land
            call get_indices2d(xnew, ynew, x0, y0, dx, dy, ig, jg)
            if (depdata(ig, jg) .lt. 0.0d0) then
              !---------------------------------------------
              ! Trying a different approach
              DBG, "Bouncing again"
              debug(seamask(ig_prev, jg_prev))
              debug(seamask(ig, jg))
              if (ig .ne. ig_prev) then
                ! Particle has moved to the left or right box
                cartx_new = cartx - 0.5 * (pvelu + pvelunew) * dt
              end if
              if (jg .ne. jg_prev) then
                ! Particle has moved to the up or down
                carty_new = carty - 0.5 * (pvelv + pvelvnew) * dt
              end if
              call xy2lonlat(cartx_new, carty_new, x0, y0, xnew, ynew)
              debug(xnow); debug(ynow); debug(znow)
              debug(xnew); debug(ynew); debug(znew)
            end if
            !---------------------------------------------
            ! Diffuse
            if (do_diffusion) call diffuse(cartx_new, carty_new, znew)
            !---------------------------------------------
            ! New coordinates and indices
            call xy2lonlat(cartx_new, carty_new, x0, y0, xnew, ynew)
            call get_indices2d(xnew, ynew, x0, y0, dx, dy, ig, jg)
            if (ig .gt. nx) xnew = lons(nx) ! TODO: These checks should go somewhere else
            if (jg .gt. ny) ynew = lats(ny) !       Is this even legal?
            if (run_3d) call particles(ipart)%check_depth(zaxdatanew_interp, ig, jg, znew, .true.)
            DBG, "Old positions 2: ", xnow, ynow, znow
            DBG, "New positions out: ", xnew, ynew, znew
          else
            !---------------------------------------------
            ! Don't move if the particle is alive, but not active
            xnew = xnow
            ynew = ynow
            znew = znow
            pvelu = 0.
            pvelv = 0.
            pvelw = 0.
            pvelunew = 0.
            pvelvnew = 0.
            pvelwnew = 0.
          end if
          !---------------------------------------------
          ! Update particles
          call particles(ipart)%update
          !---------------------------------------------
          ! Check if particle has beached or reached the boundary
          call particles(ipart)%check_beached_bdy
          if (particles(ipart)%age == max_age) call write_data_snapshot(particles(ipart), ipart)
          !---------------------------------------------
          ! Check if the particle is still alive
          if (max_age > 0) call particles(ipart)%check_age(max_age)
        end do ! End of particle loop
        DBG, "End of particle loop"
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
        itime_interp = itime_interp + 1

        if (interp_first) then
          udata_interp = udatanew_interp
          vdata_interp = vdatanew_interp
          if (run_3d) then
            wdata_interp = wdatanew_interp
            zaxdata_interp = zaxdatanew_interp
          end if
        end if

      end do ! End of interpolated time loop
      DBG, "End of interpolated time loop"

      !---------------------------------------------
      ! Set the data at the next time step as the next ... current ... data
      ! This way two fields are only read once
      ! There may be problems when doing 3D simulations (field data as derived type?)
      if (read_first) then
        udata = udatanew
        vdata = vdatanew

        if (run_3d) then
          wdata = wdatanew
          zaxdata = zaxdatanew
        end if

        if (has_viscosity) then
          visc = viscnew
        end if

        select case (has_density)
        case (DENSITY)
          density = densitynew
        case (TEMP_SALT)
          temp = tempnew
          salt = saltnew
        end select

      end if

    end do ! End of nc(t) -> nc(t+1) time loop
    DBG, "End of time loop"

    call print_loop_end(itime)

    dbgtail(loop)
  end subroutine loop
  !===========================================
  subroutine print_loop_start

    call date_and_time(date=d, time=t)
    FMT1, "======== Starting time loop ========"

    return
  end subroutine print_loop_start
  !===========================================
  subroutine print_loop_end(itime)
    integer, intent(in) :: itime

    call date_and_time(date=d, time=t)
    FMT2, LINE
    FMT2, "Finished all time steps (", itime, ")"
    FMT2, LINE

    return
  end subroutine print_loop_end

  !===========================================
  subroutine get_current_speed(xloc, yloc, zloc, uarr, varr, warr, zaxarr, uspeedout, vspeedout, wspeedout)
    !---------------------------------------------
    ! idk if this belongs here
    !---------------------------------------------

    real(rk), intent(in)            :: xloc, yloc
    real(rk), intent(in), optional :: zloc
    real(rk), intent(in)            :: uarr(nx, ny, nlevels), varr(nx, ny, nlevels)
    real(rk), intent(in), optional :: warr(nx, ny, nlevels), zaxarr(nx, ny, nlevels)
    real(rk), intent(out)           :: uspeedout, vspeedout
    real(rk), intent(out), optional :: wspeedout

    dbghead(get_current_speed)

    call get_indices2d(xloc, yloc, x0, y0, dx, dy, ig, jg, igr, jgr)

    if ((run_3d) .and. present(zloc) .and. present(warr) .and. present(wspeedout)) then
      !call get_indices_vert(zloc, depdata(ig, jg), dz, kg, kgr)
      call get_indices_vert(zaxarr=zaxarr, zin=zloc, i=ig, j=jg, k=kg, kr=kgr, dzout=dz)
      debug(ig); debug(jg); debug(igr); debug(jgr); debug(kg); debug(kgr)
      debug(dz)
      call trilinearinterp(float(ig), float(ig + 1), float(jg), float(jg + 1), &
                           float(kg), float(kg + 1), &
                           uarr(ig, jg, kg), uarr(ig + 1, jg, kg), uarr(ig, jg + 1, kg), uarr(ig + 1, jg + 1, kg), &
                           uarr(ig, jg, kg + 1), uarr(ig + 1, jg, kg + 1), uarr(ig, jg + 1, kg + 1), uarr(ig + 1, jg + 1, kg + 1), &
                           igr, jgr, kgr, uspeedout)
      call trilinearinterp(float(ig), float(ig + 1), float(jg), float(jg + 1), &
                           float(kg), float(kg + 1), &
                           varr(ig, jg, kg), varr(ig + 1, jg, kg), varr(ig, jg + 1, kg), varr(ig + 1, jg + 1, kg), &
                           varr(ig, jg, kg + 1), varr(ig + 1, jg, kg + 1), varr(ig, jg + 1, kg + 1), varr(ig + 1, jg + 1, kg + 1), &
                           igr, jgr, kgr, vspeedout)
      call trilinearinterp(float(ig), float(ig + 1), float(jg), float(jg + 1), &
                           float(kg), float(kg + 1), &
                           warr(ig, jg, kg), warr(ig + 1, jg, kg), warr(ig, jg + 1, kg), warr(ig + 1, jg + 1, kg), &
                           warr(ig, jg, kg + 1), warr(ig + 1, jg, kg + 1), warr(ig, jg + 1, kg + 1), warr(ig + 1, jg + 1, kg + 1), &
                           igr, jgr, kgr, wspeedout)
      dbgtail(get_current_speed)
      return
    end if

    call bilinearinterp(float(ig), float(ig + 1), float(ig), float(ig + 1), &
                        float(jg), float(jg + 1), &
                        uarr(ig, jg, 1), uarr(ig + 1, jg, 1), uarr(ig, jg + 1, 1), uarr(ig + 1, jg + 1, 1), &
                        igr, jgr, uspeedout)
    call bilinearinterp(float(ig), float(ig + 1), float(ig), float(ig + 1), &
                        float(jg), float(jg + 1), &
                        varr(ig, jg, 1), varr(ig + 1, jg, 1), varr(ig, jg + 1, 1), varr(ig + 1, jg + 1, 1), &
                        igr, jgr, vspeedout)

    dbgtail(get_current_speed)
    return
  end subroutine get_current_speed
  !===========================================
  subroutine get_date_and_path_dirlist
    !---------------------------------------------
    ! Find the correct path to data from directory list.
    ! Also get the first dates in the current and next folder.
    ! Called when has_subdomains=.true.
    ! TODO: Might be a bad idea to use globals.
    !---------------------------------------------

    dbghead(get_date_and_path_dirlist)

    call find_folder(theDate, thisPath)
    call nc_get_dim(trim(thisPath)//PROC0, 'time', ncNTimes)
    call find_folder(theDate%nextDate(nc_timestep * ncNTimes), nextPath)
    if (thisPath == nextPath) then
      dateThis = init_datetime_from_netcdf(trim(thisPath)//PROC0)
      dateNext = init_datetime_from_netcdf(trim(nextPath)//PROC0, ncNTimes)
      dateNext = dateNext%nextDate(nc_timestep)
      dbgtail(get_date_and_path_dirlist)
      return
    end if
    dateThis = init_datetime_from_netcdf(trim(thisPath)//PROC0)
    dateNext = init_datetime_from_netcdf(trim(nextPath)//PROC0, 1)
    !---------------------------------------------
    ! Overwrite nextPath to be the next timestep
    call find_folder(theDate%nextDate(dt), nextPath)

    dbgtail(get_date_and_path_dirlist)
  end subroutine get_date_and_path_dirlist
  !===========================================
  subroutine get_date_and_path_filelist
    !---------------------------------------------
    ! Find the correct path to data from file list.
    ! Also get the first dates in the current and next file.
    ! Called when has_subdomains=.false.
    !---------------------------------------------

    dbghead(get_date_and_path_filelist)

    call find_file(theDate, thisPath)
    call nc_get_dim(trim(thisPath), 'time', ncNTimes)
    call find_file(theDate%nextDate(nc_timestep * ncNTimes), nextPath)
    if (thisPath == nextPath) then
      DBG, "thisPath == nextPath"
      dateThis = init_datetime_from_netcdf(trim(thisPath))
      dateNext = init_datetime_from_netcdf(trim(nextPath), ncNTimes)
      dateNext = dateNext%nextDate(nc_timestep) ! Not confusing at all :)))
      dbgtail(get_date_and_path_filelist)
      return
    end if
    DBG, "thisPath /= nextPath"
    dateThis = init_datetime_from_netcdf(trim(thisPath))
    dateNext = init_datetime_from_netcdf(trim(nextPath), 1)
    !---------------------------------------------
    ! Overwrite nextPath to be the next timestep
    call find_file(theDate%nextDate(dt), nextPath)

    dbgtail(get_date_and_path_filelist)
  end subroutine get_date_and_path_filelist
end module loop_particle
