!---------------------------------------------
!For printing
#ifdef DEBUG
#define PROGRESSINFO 1
#else
#define PROGRESSINFO 1000
#endif

#define LINE '--------------------------'
#define STAR '**************************'

#ifndef WRITESTDOUT
#define FMT1 print *, '  '
#define FMT2 print *, '      '
#define FMT3 print *, '          '
#define ERROR write(*, *) achar(27)   //'[31m  ERROR: '//achar(27)//'[0m'
#define WARNING write(*, *) achar(27) //'[33m  WARNING: '//achar(27)//'[0m'
#define DBG write(*, *) '  DBG: '
#else
#define STDOUT 6
#define FMT1 write(STDOUT, *) '  '
#define FMT2 write(STDOUT, *) '      '
#define FMT3 write(STDOUT, *) '          '
#define ERROR write(STDOUT, *) '  ERROR: '
#define WARNING write(STDOUT, *) '  WARNING: '
#define DBG write(STDOUT, *) '  DBG: '
#endif

#define var2val(name) "name = ", name
#define var2val_char(name) "name = ", trim(name)

#define dbghead(name) DBG, STAR; DBG, "*** Entering name ***"
#define dbgtail(name) DBG, "*** Leaving name ***"; DBG, STAR
#define debug(var) DBG, "[var]: ", var

#ifndef DEBUG
#undef DBG
#define DBG !just a comment...
#endif

#ifdef SAY_LESS
#undef FMT1
#undef FMT2
#undef FMT3
#undef DBG
#define FMT1 !just a comment...
#define FMT2 !just a comment...
#define FMT3 !just a comment...
#define DBG !just a comment...
#endif
!---------------------------------------------
! Allow index modification
#ifdef SNAP_TO_BOUNDS
#define AIM_INTENT 
#else
! If not allowed, some subroutines should have intent(in) instead
#define AIM_INTENT , intent(in)
#endif
!---------------------------------------------
! Some shortcuts
#define PROC0 "/"//trim(file_prefix)//"0000"//trim(file_suffix)//".nc"
! GETM usually fills with 9999., could be different in other models...
#define MISSING_VAL -99.0d0
#define ZERO 0.0d0
#define HALF 0.5d0
#define ONE 1.0d0
#define SMALL 0.0001d0
