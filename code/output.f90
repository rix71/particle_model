#include "cppdefs.h"
module output
  !----------------------------------------------------------------
  ! Module for writing the output
  !----------------------------------------------------------------
  use precdefs
  use errors
  use particle_type
  use run_params, only: runid
  use particle_vars, only: particles
  use time_vars, only: theDate, nTimes, dt
  use field_vars, only: run_3d
  use netcdf
  use nc_manager, only: nc_write, nc_initialise, &
                        nc_add_dimension, nc_add_variable, &
                        nc_add_attr, nc_check, &
                        FILLVALUE_BIG
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: outputstep, snap_interval, init_output, open_beach_bdy_files, &
            close_beach_bdy_files, write_data, &
            write_beached, write_boundary, write_data_only_active, &
            write_all_particles, write_active_particles, &
            write_data_snapshot, outDir, write_snapshot
  !---------------------------------------------
  integer                     :: outputstep
  character(len=512)          :: outDir
  integer                     :: nc_t_dimid, nc_p_dimid
  character(len=512)          :: nc_fileout_all, nc_fileout_active, nc_fileout_snap
  logical                     :: write_all_particles, write_active_particles, write_snapshot
  real(rk)                    :: snap_interval
  namelist /output_vars/ outDir, outputstep, snap_interval, write_all_particles, write_active_particles, write_snapshot
  !---------------------------------------------
  integer                     :: ierr
  !===================================================
contains
  !===========================================
  subroutine init_output

    logical :: dirExists

    FMT1, "======== Init output ========"

    open (NMLFILE, file=NMLFILENAME, action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error("init_output", "Failed to open "//NMLFILENAME, ierr)
    read (NMLFILE, nml=output_vars)
    close (NMLFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error("init_output", "Failed to close "//NMLFILENAME, ierr)

    inquire (file=trim(outDir), exist=dirExists)
    if (.not. dirExists) then
#ifndef NOSYSCALLS
      FMT2, "Making directory "//trim(outDir)
      call system('mkdir -p '//trim(outDir))
#else
      call throw_error("init_output", "Out dir ("//trim(outDir)//") does not exist!")
#endif
    end if

    FMT2, var2val(write_all_particles)
    FMT2, var2val(write_active_particles)
    FMT2, var2val(write_snapshot)
    FMT2, var2val(snap_interval)
    FMT2, "Writing output every ", outputstep, " timesteps, or ", (outputstep * dt) / 3600., "hours"
    FMT2, "Saving ", nTimes / outputstep, " timesteps"

    if (write_all_particles) then
      nc_fileout_all = trim(outDir)//'/'//trim(runid)//'.all.nc'
      call init_nc_output(nc_fileout_all)
    end if

    if (write_active_particles) then
      nc_fileout_active = trim(outDir)//'/'//trim(runid)//'.active.nc'
      call init_nc_output(nc_fileout_active)
    end if

    if (write_snapshot) then
      nc_fileout_snap = trim(outDir)//'/'//trim(runid)//'.snap.nc'
      call nc_initialise(nc_fileout_snap)
      call nc_add_dimension(nc_fileout_snap, "particle", nc_p_dimid)
      call nc_add_variable(nc_fileout_snap, "time", "float", 1, [nc_p_dimid])
      call nc_add_attr(nc_fileout_snap, "time", "units", "seconds since 1900-01-01 00:00:00")

      call nc_add_variable(nc_fileout_snap, "x", "float", 1, [nc_p_dimid], FILLVALUE_BIG)
      call nc_add_attr(nc_fileout_snap, "x", "units", "degrees east")

      call nc_add_variable(nc_fileout_snap, "y", "float", 1, [nc_p_dimid], FILLVALUE_BIG)
      call nc_add_attr(nc_fileout_snap, "y", "units", "degrees north")

      call nc_add_variable(nc_fileout_snap, "z", "float", 1, [nc_p_dimid], FILLVALUE_BIG)
      call nc_add_attr(nc_fileout_snap, "z", "units", "m")
      call nc_add_attr(nc_fileout_snap, "z", "name", "depth")

      call nc_add_variable(nc_fileout_snap, "vx", "float", 1, [nc_p_dimid], FILLVALUE_BIG)
      call nc_add_attr(nc_fileout_snap, "vx", "units", "m/s")
      call nc_add_attr(nc_fileout_snap, "vx", "name", "eastward velocity")

      call nc_add_variable(nc_fileout_snap, "vy", "float", 1, [nc_p_dimid], FILLVALUE_BIG)
      call nc_add_attr(nc_fileout_snap, "vy", "units", "m/s")
      call nc_add_attr(nc_fileout_snap, "vy", "name", "northward velocity")

      call nc_add_variable(nc_fileout_snap, "vz", "float", 1, [nc_p_dimid], FILLVALUE_BIG)
      call nc_add_attr(nc_fileout_snap, "vz", "units", "m/s")
      call nc_add_attr(nc_fileout_snap, "vz", "name", "vertical velocity")

      call nc_add_variable(nc_fileout_snap, "age", "float", 1, [nc_p_dimid], FILLVALUE_BIG)
      call nc_add_attr(nc_fileout_snap, "age", "units", "s")

      call nc_add_variable(nc_fileout_snap, "id", "float", 1, [nc_p_dimid], FILLVALUE_BIG)
      call nc_add_attr(nc_fileout_snap, "id", "name", "particle id")

      call nc_add_variable(nc_fileout_snap, "trajectory", "float", 1, [nc_p_dimid], FILLVALUE_BIG)
      call nc_add_attr(nc_fileout_snap, "trajectory", "units", "m")
      call nc_add_attr(nc_fileout_snap, "trajectory", "name", "distance travelled")

      call nc_add_variable(nc_fileout_snap, "particle_num", "int", 1, [nc_p_dimid])
    end if

    return
  end subroutine init_output
  !===========================================
  subroutine init_nc_output(file_name)

    character(len=512), intent(in) :: file_name

    call nc_initialise(file_name)
    call nc_add_dimension(file_name, "particle", nc_p_dimid)
    call nc_add_dimension(file_name, "time", nc_t_dimid)

    call nc_add_variable(file_name, "time", "float", 1, [nc_t_dimid])
    call nc_add_attr(file_name, "time", "units", "seconds since 1900-01-01 00:00:00")

    call nc_add_variable(file_name, "x", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "x", "units", "degrees east")

    call nc_add_variable(file_name, "y", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "y", "units", "degrees north")

    call nc_add_variable(file_name, "z", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "z", "units", "m")
    call nc_add_attr(file_name, "z", "name", "depth")

    call nc_add_variable(file_name, "vx", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "vx", "units", "m/s")
    call nc_add_attr(file_name, "vx", "name", "eastward velocity")

    call nc_add_variable(file_name, "vy", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "vy", "units", "m/s")
    call nc_add_attr(file_name, "vy", "name", "northward velocity")

    call nc_add_variable(file_name, "vz", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "vz", "units", "m/s")
    call nc_add_attr(file_name, "vz", "name", "vertical velocity")

    call nc_add_variable(file_name, "age", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "age", "units", "s")

    call nc_add_variable(file_name, "id", "float", 1, [nc_p_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "id", "name", "particle id")

    call nc_add_variable(file_name, "trajectory", "float", 2, [nc_p_dimid, nc_t_dimid], FILLVALUE_BIG)
    call nc_add_attr(file_name, "trajectory", "units", "m")
    call nc_add_attr(file_name, "trajectory", "name", "distance travelled")

  end subroutine init_nc_output
  !===========================================
  subroutine open_beach_bdy_files

    FMT2, "Opening beached.dat..."
    open (BCHFILE, file=trim(outDir)//'/'//'beached.dat', iostat=ierr)
    if (ierr .ne. 0) call throw_error("open_beach_bdy_files", "Failed to open beached.dat", ierr)
    FMT2, "Opening boundary.dat..."
    open (BDYFILE, file=trim(outDir)//'/'//'boundary.dat', iostat=ierr)
    if (ierr .ne. 0) call throw_error("open_beach_bdy_files", "Failed to open boundary.dat", ierr)

    return
  end subroutine open_beach_bdy_files
  !===========================================
  subroutine close_beach_bdy_files

    FMT2, "Closing beached.dat..."
    close (BCHFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error("close_beach_bdy_files", "Failed to close beached.dat", ierr)
    FMT2, "Closing boundary.dat..."
    close (BDYFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error("close_beach_bdy_files", "Failed to close boundary.dat", ierr)

    return
  end subroutine close_beach_bdy_files
  !===========================================
  subroutine write_data(nwrite)
    !---------------------------------------------
    ! Write the output
    ! TODO: selection for output
    !---------------------------------------------

    integer, intent(in) :: nwrite
    integer             :: ipart, ncid, varid
    integer, save       :: nc_itime_out = 0
    real(rk)            :: var1d(1), var2d(1, 1), dateval(1)

    call theDate%print_short_date
    FMT2, "Saving data... ", nwrite, " particles"

    nc_itime_out = nc_itime_out + 1
    call nc_check(trim(nc_fileout_all), nf90_open(trim(nc_fileout_all), nf90_write, ncid), "write_data :: open")
    call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "time", varid), "write_data :: inq varid")
    dateval = theDate%date2num()
    debug(dateval)
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, dateval, start=[nc_itime_out], count=[1]), &
                  "write_data :: put var")
    do ipart = 1, nwrite
      call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "x", varid), "write_data :: inq varid")
      var2d = particles(ipart)%xPos
      call nc_check(trim(nc_fileout_all), &
                    nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                    "write_data :: put var")

      call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "y", varid), "write_data :: inq varid")
      var2d = particles(ipart)%yPos
      call nc_check(trim(nc_fileout_all), &
                    nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                    "write_data :: put var")

      call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "z", varid), "write_data :: inq varid")
      var2d = particles(ipart)%zPos
      call nc_check(trim(nc_fileout_all), &
                    nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                    "write_data :: put var")

      call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "vx", varid), "write_data :: inq varid")
      var2d = particles(ipart)%u
      call nc_check(trim(nc_fileout_all), &
                    nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                    "write_data :: put var")

      call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "vy", varid), "write_data :: inq varid")
      var2d = particles(ipart)%v
      call nc_check(trim(nc_fileout_all), &
                    nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                    "write_data :: put var")

      call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "vz", varid), "write_data :: inq varid")
      var2d = particles(ipart)%w
      call nc_check(trim(nc_fileout_all), &
                    nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                    "write_data :: put var")

      call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "age", varid), "write_data :: inq varid")
      var2d = particles(ipart)%age
      call nc_check(trim(nc_fileout_all), &
                    nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                    "write_data :: put var")

      call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "id", varid), "write_data :: inq varid")
      var1d = particles(ipart)%originNum
      call nc_check(trim(nc_fileout_all), &
                    nf90_put_var(ncid, varid, var1d, start=[ipart], count=[1]), &
                    "write_data :: put var")

      call nc_check(trim(nc_fileout_all), nf90_inq_varid(ncid, "trajectory", varid), "write_data :: inq varid")
      var2d = particles(ipart)%trajLen
      call nc_check(trim(nc_fileout_all), &
                    nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                    "write_data :: put var")
    end do
    call nc_check(trim(nc_fileout_all), nf90_close(ncid), "write :: close")

    return
  end subroutine write_data
  !===========================================
  subroutine write_data_only_active(nwrite)

    integer, intent(in) :: nwrite
    integer             :: ipart, ncid, varid
    integer, save       :: nc_itime_out = 0
    real(rk)            :: var1d(1), var2d(1, 1), dateval(1)

    nc_itime_out = nc_itime_out + 1
    call nc_check(trim(nc_fileout_active), nf90_open(trim(nc_fileout_active), nf90_write, ncid), "write_data_active :: open")
    call nc_check(trim(nc_fileout_active), nf90_inq_varid(ncid, "time", varid), "write_data_active :: inq varid")
    dateval = theDate%date2num()
    call nc_check(trim(nc_fileout_all), nf90_put_var(ncid, varid, dateval, start=[nc_itime_out], count=[1]), &
                  "write_data_active :: put var")

    do ipart = 1, nwrite
      if (particles(ipart)%isActive) then
        call nc_check(trim(nc_fileout_active), nf90_inq_varid(ncid, "x", varid), "write_data_active :: inq varid")
        var2d = particles(ipart)%xPos
        call nc_check(trim(nc_fileout_active), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data_active :: put var")

        call nc_check(trim(nc_fileout_active), nf90_inq_varid(ncid, "y", varid), "write_data_active :: inq varid")
        var2d = particles(ipart)%yPos
        call nc_check(trim(nc_fileout_active), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data_active :: put var")

        call nc_check(trim(nc_fileout_active), nf90_inq_varid(ncid, "z", varid), "write_data_active :: inq varid")
        var2d = particles(ipart)%zPos
        call nc_check(trim(nc_fileout_active), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data_active :: put var")

        call nc_check(trim(nc_fileout_active), nf90_inq_varid(ncid, "vx", varid), "write_data_active :: inq varid")
        var2d = particles(ipart)%u
        call nc_check(trim(nc_fileout_active), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data_active :: put var")

        call nc_check(trim(nc_fileout_active), nf90_inq_varid(ncid, "vy", varid), "write_data_active :: inq varid")
        var2d = particles(ipart)%v
        call nc_check(trim(nc_fileout_active), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data_active :: put var")

        call nc_check(trim(nc_fileout_active), nf90_inq_varid(ncid, "vz", varid), "write_data_active :: inq varid")
        var2d = particles(ipart)%w
        call nc_check(trim(nc_fileout_active), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data_active :: put var")

        call nc_check(trim(nc_fileout_active), nf90_inq_varid(ncid, "age", varid), "write_data_active :: inq varid")
        var2d = particles(ipart)%age
        call nc_check(trim(nc_fileout_active), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data_active :: put var")

        call nc_check(trim(nc_fileout_active), nf90_inq_varid(ncid, "id", varid), "write_data_active :: inq varid")
        var1d = particles(ipart)%originNum
        call nc_check(trim(nc_fileout_active), &
                      nf90_put_var(ncid, varid, var1d, start=[ipart], count=[1]), &
                      "write_data_active :: put var")

        call nc_check(trim(nc_fileout_active), nf90_inq_varid(ncid, "trajectory", varid), "write_data_active :: inq varid")
        var2d = particles(ipart)%trajLen
        call nc_check(trim(nc_fileout_active), &
                      nf90_put_var(ncid, varid, var2d, start=[ipart, nc_itime_out], count=[1, 1]), &
                      "write_data_active :: put var")
      end if
    end do

    call nc_check(trim(nc_fileout_active), nf90_close(ncid), "write_data_active :: close")

  end subroutine write_data_only_active
  !===========================================
  subroutine write_data_snapshot(p, particle_num)

    class(particle), intent(in) :: p
    integer, intent(in) :: particle_num
    integer             :: ncid, varid
    integer, save       :: nc_itime_out = 0
    real(rk)            :: var1d(1), dateval(1)

    nc_itime_out = nc_itime_out + 1
    call nc_check(trim(nc_fileout_snap), nf90_open(trim(nc_fileout_snap), nf90_write, ncid), "write_data_snap :: open")
    call nc_check(trim(nc_fileout_snap), nf90_inq_varid(ncid, "time", varid), "write_data_snap :: inq varid")
    dateval = theDate%date2num()

    call nc_check(trim(nc_fileout_snap), nf90_put_var(ncid, varid, dateval, start=[nc_itime_out], count=[1]), &
                  "write_data_snap :: put var")
    call nc_check(trim(nc_fileout_snap), nf90_inq_varid(ncid, "x", varid), "write_data_snap :: inq varid")
    var1d = p%xPos
    call nc_check(trim(nc_fileout_snap), &
                  nf90_put_var(ncid, varid, var1d, start=[nc_itime_out], count=[1]), &
                  "write_data_snap :: put var")

    call nc_check(trim(nc_fileout_snap), nf90_inq_varid(ncid, "y", varid), "write_data_snap :: inq varid")
    var1d = p%yPos
    call nc_check(trim(nc_fileout_snap), &
                  nf90_put_var(ncid, varid, var1d, start=[nc_itime_out], count=[1]), &
                  "write_data_snap :: put var")

    call nc_check(trim(nc_fileout_snap), nf90_inq_varid(ncid, "z", varid), "write_data_snap :: inq varid")
    var1d = p%zPos
    call nc_check(trim(nc_fileout_snap), &
                  nf90_put_var(ncid, varid, var1d, start=[nc_itime_out], count=[1]), &
                  "write_data_snap :: put var")

    call nc_check(trim(nc_fileout_snap), nf90_inq_varid(ncid, "vx", varid), "write_data_active :: inq varid")
    var1d = p%u
    call nc_check(trim(nc_fileout_snap), &
                  nf90_put_var(ncid, varid, var1d, start=[nc_itime_out], count=[1]), &
                  "write_data_snap :: put var")

    call nc_check(trim(nc_fileout_snap), nf90_inq_varid(ncid, "vy", varid), "write_data_snap :: inq varid")
    var1d = p%v
    call nc_check(trim(nc_fileout_snap), &
                  nf90_put_var(ncid, varid, var1d, start=[nc_itime_out], count=[1]), &
                  "write_data_snap :: put var")

    call nc_check(trim(nc_fileout_snap), nf90_inq_varid(ncid, "vz", varid), "write_data_snap :: inq varid")
    var1d = p%w
    call nc_check(trim(nc_fileout_snap), &
                  nf90_put_var(ncid, varid, var1d, start=[nc_itime_out], count=[1]), &
                  "write_data_snap :: put var")

    call nc_check(trim(nc_fileout_snap), nf90_inq_varid(ncid, "age", varid), "write_data_snap :: inq varid")
    var1d = p%age
    call nc_check(trim(nc_fileout_snap), &
                  nf90_put_var(ncid, varid, var1d, start=[nc_itime_out], count=[1]), &
                  "write_data_snap :: put var")

    call nc_check(trim(nc_fileout_snap), nf90_inq_varid(ncid, "id", varid), "write_data_snap :: inq varid")
    var1d = p%originNum
    call nc_check(trim(nc_fileout_snap), &
                  nf90_put_var(ncid, varid, var1d, start=[nc_itime_out], count=[1]), &
                  "write_data_snap :: put var")

    call nc_check(trim(nc_fileout_snap), nf90_inq_varid(ncid, "trajectory", varid), "write_data_snap :: inq varid")
    var1d = p%trajLen
    call nc_check(trim(nc_fileout_snap), &
                  nf90_put_var(ncid, varid, var1d, start=[nc_itime_out], count=[1]), &
                  "write_data_snap :: put var")

    call nc_check(trim(nc_fileout_snap), nf90_inq_varid(ncid, "particle_num", varid), "write_data_snap :: inq varid")
    var1d = particle_num
    call nc_check(trim(nc_fileout_snap), &
                  nf90_put_var(ncid, varid, var1d, start=[nc_itime_out], count=[1]), &
                  "write_data_snap :: put var")

    call nc_check(trim(nc_fileout_snap), nf90_close(ncid), "write_data_snap :: close")

  end subroutine write_data_snapshot
  !===========================================
  subroutine write_beached(pin)

    type(particle), intent(in) :: pin

    if (pin%state .eq. 1.) then
      write (BCHFILE, *) pin%xPos, pin%yPos, pin%originNum, pin%age, pin%trajLen
    else
      call throw_warning("write_beached", "This particle should not be beached")
    end if

    return
  end subroutine write_beached
  !===========================================
  subroutine write_boundary(pin)

    type(particle), intent(in) :: pin

    if (pin%state .eq. 2.) then
      write (BDYFILE, *) pin%xPos, pin%yPos, pin%originNum, pin%age, pin%trajLen
    else
      call throw_warning("write_boundary", "This particle should not be on boundary")
    end if

    return
  end subroutine write_boundary

end module output
