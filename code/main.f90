#include "cppdefs.h"
#include "output.h"
program main
  use mod_errors
  use mod_precdefs
#ifdef WRITESTDOUT
  use run_params, only: runid
#endif
  use run_params, only: dry_run, nmlfilename
! TODO: General biofouling flag
#if (defined(BIOFOULING_KOOI) || defined(BIOFOULING_TSIARAS) || defined(BIOFOULING_SIMPLE))
  use run_params, only: biofouling_nmlfilename
#endif
  use mod_initialise, only: init_run, init_model
  use mod_loop, only: loop
  ! use mod_output, only: init_output
#ifdef POSTPROCESS
  use mod_postprocessing, only: postprocess
#endif
  implicit none
  !===================================================
  call command_line
  call init_run
#ifdef WRITESTDOUT
  open (STDOUT, file=trim(runid)//".stdout")
#endif
  call init_model
  ! call init_output
  !---------------------------------------------
  if (.not. dry_run) then
    call loop
#ifdef POSTPROCESS
    call postprocess
#endif
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
      case ('--variables')
        call print_variables()
        stop
      case ('-nml', '--namelist')
        call getarg(i + 1, nmlfilename)
        i = i + 1
#if (defined(BIOFOULING_KOOI) || defined(BIOFOULING_TSIARAS) || defined(BIOFOULING_SIMPLE))
      case ('-bnml', '--biofouling-namelist')
        call getarg(i + 1, biofouling_nmlfilename)
        i = i + 1
#endif
      case default
        ERROR, "Command line argument not recognised: "//trim(arg)
        call print_help()
        stop
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
    FMT1, "Command line options:"
    FMT2, "-h, --help                               print this help message and exit"
    FMT2, "-c, --compile                            print compliation options and exit"
    FMT2, "--variables                              print output variables and exit"
    FMT2, "-nml, --namelist <filename>              use <filename> as namelist (default: input.inp)"
#if (defined(BIOFOULING_KOOI) || defined(BIOFOULING_TSIARAS) || defined(BIOFOULING_SIMPLE))
    FMT2, "-bnml, --biofouling-namelist <filename>  use <filename> as biofouling namelist (default: biofouling.inp)"
#endif
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
#ifdef NO_ADVECT_VERTICAL
    FMT2, "-NO_ADVECT_VERTICAL"
#endif
#ifdef NO_DIFFUSE_VERTICAL
    FMT2, "-NO_DIFFUSE_VERTICAL"
#endif
#ifdef SMAGORINSKY_FULL_FIELD
    FMT2, "-SMAGORINSKY_FULL_FIELD"
#endif
#ifdef SMAGORINSKY_INTERP_UV
    FMT2, "-SMAGORINSKY_INTERP_UV"
#endif
#ifdef POSTPROCESS
    FMT2, "-POSTPROCESS"
#endif
#ifdef BIOFOULING_KOOI
    FMT2, "-BIOFOULING_KOOI"
#endif
#ifdef BIOFOULING_SIMPLE
    FMT2, "-BIOFOULING_SIMPLE"
#endif
#ifdef BIOFOULING_TSIARAS
    FMT2, "-BIOFOULING_TSIARAS"
#endif

    return
  end subroutine print_compile_info
  !===========================================
  subroutine print_variables
    FMT1, ""
    FMT1, "Output variables:"
    FMT1, LINE
    FMT1, "time:           time"
    FMT1, "lon:            longitude"
    FMT1, "lat:            latitude"
    FMT1, "depth:          depth"
#ifdef OUT_ID
    FMT1, "id:             particle (source) id"
#endif
#ifdef OUT_VELOCITY
    FMT1, "vx:             x-velocity"
    FMT1, "vy:             y-velocity"
    FMT1, "vz:             z-velocity"
#endif
#ifdef OUT_SETTLING_VELOCITY
    FMT1, "vs:             settling velocity"
#endif
#ifdef OUT_DENSITY
    FMT1, "rho:            density"
#endif
#ifdef OUT_DENSITY_PLASTIC
    FMT1, "rho_plastic:    density of plastic particle"
#endif
#ifdef OUT_DELTA_RHO
    FMT1, "delta_rho:      density difference between particle and fluid"
#endif
#ifdef OUT_RADIUS
    FMT1, "radius:         radius"
#endif
#ifdef OUT_RADIUS_PLASTIC
    FMT1, "radius_plastic: radius of plastic particle"
#endif
#ifdef OUT_AGE
    FMT1, "age:            age of particle"
#endif
#ifdef OUT_TRAJECTORY
    FMT1, "trajectory:     total distance travelled by particle"
#endif
#ifdef OUT_TIME_ON_BEACH
    FMT1, "time_on_beach:  time particle has been on beach"
#endif
#ifdef OUT_BEACHING_TIME
    FMT1, "beaching_time:  maximum time particle can be on beach before it is removed"
#endif
#ifdef OUT_KIN_VISCOSITY
    FMT1, "kin_viscosity:  kinematic viscosity"
#endif
#ifdef OUT_FRICTION_VELOCITY
    FMT1, "u_star:         friction velocity"
#endif
#ifdef OUT_H_BIOFILM
    FMT1, "h_biofilm:      thickness of biofilm"
#endif
#ifdef OUT_GROWTH_BIOFILM
    FMT1, "growth_biofilm: growth rate of biofilm"
#endif
#ifdef OUT_STATE 
    FMT1, "state:          state of particle"
#endif
  end subroutine print_variables
end program main
