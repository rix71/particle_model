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
        call print_help()
        stop
      case ('-c', '--compile')
        call print_compile_info()
        stop
      case default
        call print_help()
        call throw_error("main :: command_line", "Command line argument not recognised: "//trim(arg))
      end select
      i = i + 1
    end do

    return
  end subroutine command_line
  !===========================================
#ifdef SAY_LESS
#define FMT1 print *, '  '
#define FMT2 print *, '      '
#define FMT3 print *, '          '
#endif  
  subroutine print_help()
    character(len=LEN_CHAR_L) :: cmd

    call get_command_argument(0, cmd)

    FMT1, ""
    FMT1, "Usage: ", trim(cmd), " [OPTIONS]"
    FMT1, ""
    FMT1, "Without any options, the program will continue execution."
    FMT1, "All options for execution (and dry run) are set in the namelist (name in cppdefs.h)."
    FMT1, ""
    FMT1, "Command line options:"
    FMT2, "-h, --help      print this help message and exit"
    FMT2, "-c, --compile   print compliation options and exit"
    FMT1, ""

    return
  end subroutine print_help
  !===========================================
  subroutine print_compile_info()
    use iso_fortran_env

    FMT1, ""
    FMT1, "Compiled: "//__DATE__//" "//__TIME__
    FMT1, "Compiler: "//compiler_version()
    FMT1, "Compiler options: "//compiler_options()
    FMT1, ""
    FMT1, "Used flags:"
#ifdef WRITESTDOUT
    FMT2, "-WRITESTDOUT"
#endif
#ifdef USE_OMP
    FMT2, "-USE_OMP"
#endif
#ifdef DEBUG
    FMT2, "-DEBUG"
#endif
#ifdef SAY_LESS
    FMT2, "-SAY_LESS"
#endif
#ifdef NOSYSCALLS
    FMT2, "-NOSYSCALLS"
#endif
#ifdef SNAP_TO_BOUNDS
    FMT2, "-SNAP_TO_BOUNDS"
#endif
#ifdef PARTICLE_BOUNCE
    FMT2, "-PARTICLE_BOUNCE"
#endif
#ifdef PARTICLE_REDIRECT
    FMT2, "-PARTICLE_REDIRECT"
#endif
#ifdef PARTICLE_BEACH_IMMEDIATELY
    FMT2, "-PARTICLE_BEACH_IMMEDIATELY"
#endif
#ifdef PARTICLE_SNAP_SEALVL
    FMT2, "-PARTICLE_SNAP_SEALVL"
#endif
#ifdef ADVECT_VERTICAL
    FMT2, "-ADVECT_VERTICAL"
#endif
#ifdef DIFFUSE_VERTICAL
    FMT2, "-DIFFUSE_VERTICAL"
#endif
#ifdef SMAGORINSKY_FULL_FIELD
    FMT2, "-SMAGORINSKY_FULL_FIELD"
#endif
#ifdef SMAGORINSKY_INTERP_UV
    FMT2, "-SMAGORINSKY_INTERP_UV"
#endif

    return
  end subroutine print_compile_info
end program main
