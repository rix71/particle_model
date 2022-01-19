#include "cppdefs.h"
program main
  use precdefs
#ifdef WRITESTDOUT
  use run_params, only: runid
#endif
  use run_params, only: dry_run
  use initialise, only: init_run, init_model
  use loop_particle, only: loop
  use output, only: init_output, open_beach_bdy_files, close_beach_bdy_files
  use postprocessing, only: postprocess
  implicit none
  !===================================================
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

end program main
