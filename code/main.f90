#include "cppdefs.h"
program main
  use mod_errors
  use mod_precdefs
#ifdef WRITESTDOUT
  use run_params, only: runid
#endif
  use run_params, only: dry_run
  use mod_initialise, only: init_run, init_model
  use mod_loop, only: loop
  use mod_output, only: init_output, open_beach_bdy_files, close_beach_bdy_files
  use mod_postprocessing, only: postprocess
  implicit none
  !===================================================
  call command_line()
  call init_run
#ifdef WRITESTDOUT
  open (STDOUT, file=trim(runid)//".stdout")
#endif
  call init_model
  call init_output
  !---------------------------------------------
  if (.not. dry_run) then
    call loop
    call postprocess
  else
    FMT1, LINE; FMT1, "Will not loop!"; FMT1, LINE
  end if
  !---------------------------------------------
  FMT1, "FINISHED"
#ifdef WRITESTDOUT
  close (STDOUT)
#endif
  !===================================================
contains
  !===========================================
  subroutine command_line()
    integer :: i
    character(len=256) :: arg

    i = 1
    do while (i <= iargc())
      call getarg(i, arg)
      select case (arg)
      case ('-h', '--help')
        call throw_error("main :: command_line", "Sorry, no help yet")
      case ('-c')
        call print_compile_info()
        stop
      case default
        call throw_error("main :: command_line", "Command line argument not recognised: "//trim(arg))
      end select
      i = i + 1
    end do

    return
  end subroutine command_line
  !===========================================
  subroutine print_compile_info()

    FMT1, "Compiled: "//__DATE__//" "//__TIME__

    return
  end subroutine print_compile_info
end program main
