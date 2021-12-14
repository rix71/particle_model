#include "cppdefs.h"
module loop_particle
  !----------------------------------------------------------------
  ! Main loop
  !----------------------------------------------------------------
  use precdefs
  use errors
  use interp
  use params, only: do_diffusion, do_velocity
  use physics, only: diffuse, velocity, vertical_velocity
  use loop_vars, only: ipart, itime, dateThis, dateNext, thisPath, nextPath, &
                       xnow, xnew, ynow, ynew, znow, znew, cartx, cartx_new, &
                       carty, carty_new, ig, jg, kg, igr, jgr, kgr, &
                       nc_itime, nc_itime_next, ncNTimes, pvelu, pvelunew, &
                       pvelv, pvelvnew, pvelw, pvelwnew
  use particle_type, only: particle
  use particle_vars, only: nInitParticles, runparts, inputstep, initCoords, particles
  use time_vars, only: theDate, run_end_dt, dt, nc_timestep, nTimes
  use domain, only: lonlat2xy, xy2lonlat
  use domain_vars, only: x0, y0, dx, dy, dz, dx_m, dy_m, nx, ny, lons, lats, depdata
  use field_vars, only: nlevels, has_subdomains, run_3d, read_first, &
                        uspeed, vspeed, wspeed, &
                        uspeednew, vspeednew, wspeednew, &
                        udata, vdata, wdata, &
                        udatanew, vdatanew, wdatanew, &
                        zaxdata, zaxdatanew, rho_sea, viscosity, &
                        file_prefix, file_suffix
  use fields, only: find_folder, find_file, get_indices2d, &
                    get_indices_vert, read_fields_full_domain, &
                    get_seawater_density, get_seawater_viscosity
  use nc_manager, only: nc_get_dim
  use modtime
  use output, only: outputstep, write_data, &
                    write_beached, write_boundary
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: loop
  !---------------------------------------------
  integer           :: active_particles
  character(len=8)  :: d
  character(len=10) :: t
  !===================================================
contains
  !===========================================
  subroutine loop

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
    ! Start time loop
    do while (theDate <= run_end_dt)
      if (mod(itime, PROGRESSINFO) .eq. 0) then
        call theDate%print_short_date
        call date_and_time(date=d, time=t)
        FMT2, t(1:2), ":", t(3:4), ":", t(5:10), " itime = ", itime
        if (itime .ne. 0) then
          active_particles = 0
          do ipart = 1, runparts
            if (particles(ipart)%isActive) active_particles = active_particles + 1
          end do
          FMT2, "There are ", active_particles, " active particles"
        end if
      end if
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
      ! Release particles
      ! TODO: Do this only in case of continuous release (particle_method ?)
      if (mod(itime, inputstep) .eq. 0) then
        FMT2, "Initialising particles at itime = ", itime
        do ipart = 1, nInitParticles
          particles(ipart + runparts) = particle(xPos=initCoords(ipart, 1), &
                                                 yPos=initCoords(ipart, 2), &
                                                 originNum=initCoords(ipart, 3), &
                                                 beachingtime=initCoords(ipart, 4), &
                                                 rho=initCoords(ipart, 5), &
                                                 radius=initCoords(ipart, 6))
        end do
        runparts = runparts + nInitParticles
        FMT2, runparts, "particles"
      end if
      !---------------------------------------------
      ! Get netCDF time index
      debug(thisPath)
      debug(nextPath)
      nc_itime = int(dateDiff(dateThis, theDate) / nc_timestep)
      nc_itime_next = nc_itime + 1
      debug(nc_itime)
      debug(nc_itime_next)
      debug(ncNTimes)
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
      end if
      debug(thisPath)
      debug(nextPath)
      !---------------------------------------------
      ! Read in currents
      ! Does this even make sense???
      if (runparts .gt. 0) then
        if (.not. read_first) then
          select case (run_3d)
          case (.true.)
            call read_fields_full_domain(datau=udata, datav=vdata, &
                                         dataw=wdata, datazax=zaxdata, &
                                         timeindex=nc_itime, path=thisPath)
            call read_fields_full_domain(datau=udatanew, datav=vdatanew, &
                                         dataw=wdatanew, datazax=zaxdatanew, &
                                         timeindex=nc_itime_next, path=nextPath)

            read_first = .true.
          case (.false.)
            call read_fields_full_domain(datau=udata, datav=vdata, &
                                         timeindex=nc_itime, path=thisPath)
            call read_fields_full_domain(datau=udatanew, datav=vdatanew, &
                                         timeindex=nc_itime_next, path=nextPath)
            read_first = .true.
          end select
        else
          select case (run_3d)
          case (.true.)
            call read_fields_full_domain(datau=udatanew, datav=vdatanew, &
                                         dataw=wdatanew, datazax=zaxdatanew, &
                                         timeindex=nc_itime_next, path=nextPath)
          case (.false.)
            call read_fields_full_domain(datau=udatanew, datav=vdatanew, &
                                         timeindex=nc_itime_next, path=nextPath)
          end select
        end if
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
        ! Cartesian coordinates
        call lonlat2xy(xnow, ynow, x0, y0, cartx, carty)
        !---------------------------------------------
        ! Global i and j indices
        call get_indices2d(xnow, ynow, x0, y0, dx, dy, ig, jg)
        !---------------------------------------------
        ! If particle hovers, set the depth to sealevel
        if (run_3d) then
          if (znow > zaxdata(ig, jg, nlevels)) then
            DBG, "Setting particle to sealevel"
            znow = zaxdata(ig, jg, nlevels)
          end if
        end if
        !---------------------------------------------
        ! Particle should not be on land
        if (depdata(ig, jg) .lt. 0.0d0) then
          !---------------------------------------------
          ! If depdata is -10.0 then there will be NaNs. Skip the particle.
          call throw_warning("Loop", "Particle is on land! Skipping.")
          particles(ipart)%warnings = particles(ipart)%warnings + 1
          if (particles(ipart)%warnings .ge. 3) particles(ipart)%isActive = .false.
          cycle
        end if
        !---------------------------------------------
        ! Get current speed at particle location
        select case (run_3d)
        case (.true.)
          call get_current_speed(xloc=xnow, yloc=ynow, zloc=znow, &
                                 uarr=udata, varr=vdata, warr=wdata, zaxarr=zaxdata, &
                                 uspeedout=uspeed, vspeedout=vspeed, wspeedout=wspeed)
          DBG, "Speeds 3D:"
          debug(uspeed); debug(vspeed); debug(wspeed)
        case (.false.)
          call get_current_speed(xloc=xnow, yloc=ynow, &
                                 uarr=udata, varr=vdata, &
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
          call get_seawater_density(ig, jg, kg, znow, nc_itime, rho_sea, thisPath)
          debug(rho_sea)
          !---------------------------------------------
          ! Get seawater viscosity
          call get_seawater_viscosity(ig, jg, kg, nc_itime, viscosity, thisPath)
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
          znew = znow + pvelw * dt ! TODO: Don't let the particle jump out of water
          if (znew .gt. 0.0d0) then
            DBG, "Setting znew to 0"; debug(znew)
            znew = 0.0d0 ! TODO: Which way is up???
          end if
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
          ! No point skipping now. Just keep old position
          ! TODO: Should try to decide if the particle has beached
          ! Still skip actually...
          call throw_warning("Loop", &
                             "New location of particle is on land! Keeping old position.")
          particles(ipart)%warnings = particles(ipart)%warnings + 1
          if (particles(ipart)%warnings .ge. 3) particles(ipart)%isActive = .false.
          cycle
        end if
        !---------------------------------------------
        ! Get current speed at particle location
        select case (run_3d)
        case (.true.)
          call get_current_speed(xloc=xnew, yloc=ynew, zloc=znew, &
                                 uarr=udatanew, varr=vdatanew, warr=wdatanew, zaxarr=zaxdatanew, &
                                 uspeedout=uspeednew, vspeedout=vspeednew, wspeedout=wspeednew)
          DBG, "Speeds 3D new:"
          debug(uspeednew); debug(vspeednew); debug(wspeednew)
        case (.false.)
          call get_current_speed(xloc=xnew, yloc=ynew, &
                                 uarr=udatanew, varr=vdatanew, &
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
          call get_seawater_density(ig, jg, kg, znew, nc_itime_next, rho_sea, thisPath)
          !---------------------------------------------
          ! Get seawater viscosity
          call get_seawater_viscosity(ig, jg, kg, nc_itime_next, viscosity, thisPath)
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
        ! Diffuse
        if (do_diffusion) call diffuse(cartx_new, carty_new, znew)
        !---------------------------------------------
        ! New coordinates and indices
        call xy2lonlat(cartx_new, carty_new, x0, y0, xnew, ynew)
        call get_indices2d(xnew, ynew, x0, y0, dx, dy, ig, jg)
        if (ig .gt. nx) xnew = lons(nx) ! TODO: These checks should go somewhere else
        if (jg .gt. ny) ynew = lats(ny) !       Is this even legal?
        DBG, "Old positions 2: ", xnow, ynow, znow
        DBG, "New positions out: ", xnew, ynew, znew
        !---------------------------------------------
        ! Update particles
        call particles(ipart)%update
        !---------------------------------------------
        ! Check if particle has beached or reached the boundary
        call particles(ipart)%check_beached_bdy
        if (particles(ipart)%state .eq. 1.) call write_beached(particles(ipart))
        if (particles(ipart)%state .eq. 2.) call write_boundary(particles(ipart))

      end do ! End of particle loop
      DBG, "End of particle loop"
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
      end if
      !---------------------------------------------
      ! Write output
      if ((mod(itime, outputstep) .eq. 0) .and. (runparts .gt. 0)) call write_data(itime, runparts)
      !---------------------------------------------
      ! Update time
      call theDate%update(dt)
      itime = itime + 1

    end do ! End of time loop
    DBG, "End of time loop"

    call print_loop_end

    dbgtail(loop)
  end subroutine loop
  !===========================================
  subroutine print_loop_start

    call date_and_time(date=d, time=t)
    FMT1, "======== Starting time loop ========"
    FMT2, LINE
    FMT2, "Date and time now: "
    FMT3, d(1:4), "-", d(5:6), "-", d(7:8)
    FMT3, t(1:2), ":", t(3:4), ":", t(5:10)
    FMT2, "There should be ", nTimes, "iterations"
    FMT2, "Writing output every ", outputstep, " timesteps, or ", (outputstep * dt) / 3600., "hours"
    FMT2, "This will result in ", nTimes / outputstep, " files"
    FMT2, LINE

    return
  end subroutine print_loop_start
  !===========================================
  subroutine print_loop_end

    call date_and_time(date=d, time=t)
    FMT2, LINE
    FMT2, "Finished all time steps (", itime, ")"
    FMT2, "Date and time now: "
    FMT3, d(1:4), "-", d(5:6), "-", d(7:8)
    FMT3, t(1:2), ":", t(3:4), ":", t(5:10)
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

    call find_file(theDate, thisPath)
    call nc_get_dim(trim(thisPath), 'time', ncNTimes)
    call find_file(theDate%nextDate(nc_timestep * ncNTimes), nextPath)
    if (thisPath == nextPath) then
      dateThis = init_datetime_from_netcdf(trim(thisPath))
      dateNext = init_datetime_from_netcdf(trim(nextPath), ncNTimes)
      dateNext = dateNext%nextDate(nc_timestep) ! Not confusing at all :)))
      return
    end if
    dateThis = init_datetime_from_netcdf(trim(thisPath))
    dateNext = init_datetime_from_netcdf(trim(nextPath), 1)

  end subroutine get_date_and_path_filelist
end module loop_particle
