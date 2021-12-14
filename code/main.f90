#include "cppdefs.h"
program main
  use precdefs
#ifdef WRITESTDOUT
  use run_params, only: runid
#endif
  use initialise, only: init_run, init_model
  use loop_particle, only: loop
  use output, only: init_output, open_beach_bdy_files, close_beach_bdy_files
  implicit none
  !===================================================
  call init_run
#ifdef WRITESTDOUT
  !open(STDOUT, file=STDOUTFILE)
  open (STDOUT, file=trim(runid)//".stdout")
#endif
  call init_model
  call init_output
  call open_beach_bdy_files
  !---------------------------------------------
  call loop
  !---------------------------------------------
  call close_beach_bdy_files
  FMT1, "FINISHED"
#ifdef WRITESTDOUT
  close (STDOUT)
#endif

end program main
