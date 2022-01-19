!---------------------------------------------
! For printing
#define PROGRESSINFO 1
#define LINE '--------------------------'
#define STAR '**************************'

#ifndef WRITESTDOUT 
#define FMT1 print *, '  '
#define FMT2 print *, '      '
#define FMT3 print *, '          '
#define ERROR write(*, *) achar(27)   //'[31m  ERROR: '//achar(27)//'[0m'
#define WARNING write(*, *) achar(27) //'[33m  WARNING: '//achar(27)//'[0m'
#define DBG write(*, *) '  DBG: '
#endif

#ifdef WRITESTDOUT
#define STDOUT 6
#define FMT1 write(STDOUT,*) '  '
#define FMT2 write(STDOUT,*) '      '
#define FMT3 write(STDOUT,*) '          '
#define ERROR write(STDOUT,*) '  ERROR: '
#define WARNING write(STDOUT,*) '  WARNING: '
#define DBG write(STDOUT,*) '  DBG: '
#endif

#define var2val(name) "name = ", name
#define var2val_char(name) "name = ", trim(name)

#define dbghead(name) DBG, STAR; DBG, "*** Entering name ***"
#define dbgtail(name) DBG, "*** Leaving name ***"; DBG, STAR
#define debug(var) DBG, "[var]: ", var

#ifndef DEBUG
#undef DBG
#define DBG ! just a comment...
#endif

#ifdef SAYLESS
#undef FMT1
#undef FMT2
#undef FMT3
#define FMT1 ! just a comment...
#define FMT2 ! just a comment...
#define FMT3 ! just a comment...
#endif
!---------------------------------------------
! File numbers
#define NMLFILE 10
#define COORDFILE 11
#define PROCFILE 12
#define DIRFILE 13
#define DATAOUTFILE 21
#define BCHFILE 22
#define BDYFILE 23
!---------------------------------------------
! Some shortcuts
#define PROC0 "/"//trim(file_prefix)//"0000"//trim(file_suffix)//".nc"
#define NMLFILENAME 'input.inp'
! GETM usually fills with -9999., could be different in other models...
#define MISSING_VAL -99.0d0
!---------------------------------------------
! Enumerations
!--------------
! zax_style
#define DEPTH_VALUES 1
#define LAYER_THICKNESS 2
!--------------
! has_density
#define DEFAULT_DENSITY 0
#define DENSITY 1
#define TEMP_SALT 2

