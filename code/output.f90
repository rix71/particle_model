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
  use time_vars, only: theDate
  use field_vars, only: run_3d
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: outputstep, init_output, open_beach_bdy_files, &
            close_beach_bdy_files, write_data, &
            write_beached, write_boundary
  !---------------------------------------------
  integer                     :: outputstep
  character(len=512)          :: outDir
  character(len=*), parameter :: dataDir = 'data'
  namelist /output_vars/ outDir, outputstep
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

#ifndef NOSYSCALLS
    inquire (file=trim(outDir), exist=dirExists)
    if (.not. dirExists) then
      FMT2, "Making directory "//trim(outDir)
      !call system('mkdir '//trim(outDir))
      call system('mkdir -p '//trim(outDir)//'/'//trim(dataDir))
    end if
#endif

    return
  end subroutine init_output
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
  subroutine write_data(itime, nwrite)
    !---------------------------------------------
    ! Write the output
    ! TODO: selection for output
    !---------------------------------------------

    integer, intent(in) :: itime, nwrite
    integer             :: ipart
    character(len=14)   :: outfilesuf
    character(len=512)  :: fileout

    FMT2, "Saving data... ", nwrite, " particles"

    write (outfilesuf, '(i0.14)') theDate%shortDate(.true.)
    fileout = trim(outDir)//'/'//trim(dataDir)//'/'//trim(runid)//'.'//outfilesuf//'.dat'
    open (DATAOUTFILE, file=trim(fileout), iostat=ierr)
    if (ierr .ne. 0) call throw_error("write_data", "Failed to open "//trim(fileout), ierr)
    select case (run_3d)
    case (.true.)
      write (DATAOUTFILE, *) "timestep, lon, lat, source, age, depth, trajectory"
      do ipart = 1, nwrite
        write (DATAOUTFILE, *) itime, particles(ipart)%xPos, particles(ipart)%yPos, &
          particles(ipart)%originNum, particles(ipart)%age, &
          particles(ipart)%zPos, particles(ipart)%trajLen
      end do
    case (.false.)
      write (DATAOUTFILE, *) "timestep, lon, lat, source, age, trajectory"
      do ipart = 1, nwrite
        write (DATAOUTFILE, *) itime, particles(ipart)%xPos, particles(ipart)%yPos, &
          particles(ipart)%originNum, particles(ipart)%age, &
          particles(ipart)%trajLen
      end do
    end select
    close (DATAOUTFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error("write_data", "Failed to close "//trim(fileout), ierr)

    return
  end subroutine write_data
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
