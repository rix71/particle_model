#include "cppdefs.h"
program main
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
