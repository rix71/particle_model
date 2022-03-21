#include "cppdefs.h"
module mod_datetime
  !----------------------------------------------------------------
  ! Module for handling dates
  ! Timezones are not taken into account, assumes UTC
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  use nc_manager, only: nc_get_timeunit, nc_read_time_val
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_datetime, date_diff, datetime_from_netcdf
  !---------------------------------------------
  type t_datetime
    integer           :: year = 1, month = 1, day = 1
    integer           :: hour = 0, minute = 0, second = 0
    character(len=LEN_CHAR_S) :: dtName ! Absolutely unnecessary...

  contains
    procedure :: setName
    procedure :: print_date
    procedure :: print_short_date
    procedure :: update
    procedure :: nextDate
    procedure :: isleap
    procedure :: yearday
    procedure :: date2num
    procedure :: shortDate
    !---------------------------------------------
    ! Operator overloading
    procedure, pass(this) :: date_gt
    generic :: operator(>) => date_gt
    generic :: operator(.gt.) => date_gt
    procedure, pass(this) :: date_lt
    generic :: operator(<) => date_lt
    generic :: operator(.lt.) => date_lt
    procedure, pass(this) :: date_eq
    generic :: operator(==) => date_eq
    generic :: operator(.eq.) => date_eq
    procedure, pass(this) :: date_ge
    generic :: operator(>=) => date_ge
    generic :: operator(.ge.) => date_ge
    procedure, pass(this) :: date_le
    generic :: operator(<=) => date_le
    generic :: operator(.le.) => date_le
  end type t_datetime
  !---------------------------------------------
  interface t_datetime
    procedure :: ctor_datetime
  end interface t_datetime
  !---------------------------------------------
  ! Days in month. This is pretty dumb...
  integer, dimension(12) :: daysInMonth = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
  !---------------------------------------------
  ! For conversions
  real(rk), parameter    :: d2s = 86400.
  real(rk), parameter    :: h2s = 3600.
  real(rk), parameter    :: h2d = 1./24.
  real(rk), parameter    :: min2s = 60.
  real(rk), parameter    :: min2d = 1./(24.*60.)
  real(rk), parameter    :: sec2d = 1./(24.*3600.)
  !===================================================
contains
  !===========================================
  type(t_datetime) function ctor_datetime(date_str) result(res)
    !---------------------------------------------
    ! Constructor for datetime
    ! TODO: Date validation
    !---------------------------------------------
    character(len=LEN_CHAR_S), intent(in) :: date_str
    integer                       :: year, month, day
    integer                       :: hour, minute, second
    integer                       :: len_str

    len_str = len(trim(date_str))

    if (len_str .eq. 10) then
      read (date_str, '(i4,1x,i2,1x,i2)') year, month, day
      res%year = year
      res%month = month
      res%day = day
    else if (len_str .eq. 19) then
      read (date_str, '(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)') year, month, day, &
        hour, minute, second
      res%year = year
      res%month = month
      res%day = day
      res%hour = hour
      res%minute = minute
      res%second = second
    else
      call throw_error("datetime", &
                       "Wrong date format: "//trim(date_str)// &
                       " Accepted date formats are 'yyyy-mm-dd' or 'yyyy-mm-dd HH:MM:SS'")
    end if

  end function ctor_datetime
  !===========================================
  type(t_datetime) function datetime_from_netcdf(fname, n)
    !---------------------------------------------
    ! Reads the time unit from netCDF and returns a datetime
    ! datetime instance. Time unit date if no 'n', first date
    ! in netCDF if 'n' is specified.
    !---------------------------------------------

    integer                       :: year, month, day
    integer                       :: hour, minute, second
    integer, intent(in), optional :: n
    character(len=60)             :: timeunit
    character(len=LEN_CHAR_S)             :: ncdatestr
    character(len=*), intent(in)  :: fname
    real(rk)                      :: diff

    call nc_get_timeunit(trim(fname), timeunit)
    read (timeunit, '(8x,1x,5x,i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)') year, month, day, hour, minute, second
    write (ncdatestr, "(i4,a,i2.2,a,i2.2,1x,i2.2,a,i2.2,a,i2.2)") year, "-", month, "-", day, &
      hour, ":", minute, ":", second
    datetime_from_netcdf = t_datetime(ncdatestr)

    if (present(n)) then
      diff = nc_read_time_val(trim(fname), n)
      call datetime_from_netcdf%update(diff)
    end if

  end function datetime_from_netcdf
  !===========================================
  subroutine setName(this, name)

    class(t_datetime), intent(inout) :: this
    character(len=*), intent(in)     :: name

    this%dtName = trim(name)

  end subroutine setName
  !===========================================
  subroutine print_date(this)

    class(t_datetime), intent(in) :: this

    FMT1, "============================="
    FMT1, '----', this%dtName, '----'
    FMT1, "Year: ", this%year
    FMT1, "Month: ", this%month
    FMT1, "Day: ", this%day
    FMT1, "Hour: ", this%hour
    FMT1, "Minute: ", this%minute
    FMT1, "Second: ", this%second
    FMT1, "============================="

    return
  end subroutine print_date
  !===========================================
  subroutine print_short_date(this)

    character(len=19)           :: datestr
    class(t_datetime), intent(in) :: this

    write (datestr, "(i4,a,i2.2,a,i2.2,1x,i2.2,a,i2.2,a,i2.2)") &
      this%year, "-", this%month, "-", this%day, &
      this%hour, ":", this%minute, ":", this%second
    FMT1, datestr

  end subroutine print_short_date
  !===========================================
  subroutine update(this, dt)
    !---------------------------------------------
    ! Update the date
    ! Timestep must be in seconds
    ! TODO: Backwards update
    !---------------------------------------------

    class(t_datetime), intent(inout) :: this
    real(rk), intent(in)           :: dt

    call reset_DIM
    call this%isleap
    this%second = this%second + int(dt)
    do while (this%second .ge. 60)
      this%second = this%second - 60
      this%minute = this%minute + 1
      if (this%minute .ge. 60) then
        this%minute = this%minute - 60
        this%hour = this%hour + 1
        if (this%hour .ge. 24) then
          this%hour = this%hour - 24
          this%day = this%day + 1
          if (this%day .gt. daysInMonth(this%month)) then
            this%day = this%day - daysInMonth(this%month)
            this%month = this%month + 1
            if (this%month .gt. 12) then
              this%month = this%month - 12
              this%year = this%year + 1
              call reset_DIM
              call this%isleap
            end if
          end if
        end if
      end if
    end do
    call reset_DIM

    return
  end subroutine update
  !===========================================
  type(t_datetime) function nextDate(this, dt)
    !---------------------------------------------
    ! Get a new instance of datetime at the next timestep
    ! without updating the original date
    !---------------------------------------------
    class(t_datetime), intent(in) :: this
    real(rk), intent(in)        :: dt

    nextDate = this
    call nextDate%update(dt)

  end function nextDate
  !===========================================
  subroutine isleap(this)
    !---------------------------------------------
    ! Changes the second number in daysInMonth to 29
    ! if the year is a leap year.
    ! daysInMonth then stays like this, so it must be reset
    ! TODO: Make into pure function or something so won't have to reset
    !---------------------------------------------

    class(t_datetime), intent(in) :: this

    if (mod(this%year, 4) .eq. 0) then
      daysInMonth(2) = 29
      if ((mod(this%year, 100) .eq. 0) .and. (mod(this%year, 400) .ne. 0)) then
        daysInMonth(2) = 28
      end if
    end if

    return
  end subroutine isleap
  !===========================================
  elemental logical function isLeapBool(year)
    !---------------------------------------------
    ! Check if year is a leap year
    !---------------------------------------------
    integer, intent(in) :: year

    isLeapBool = .false.
    if (mod(year, 4) .eq. 0) then
      isLeapBool = .true.
      if ((mod(year, 100) .eq. 0) .and. (mod(year, 400) .ne. 0)) then
        isLeapBool = .false.
      end if
    end if

  end function isLeapBool
  !===========================================
  integer function yearday(this)
    !---------------------------------------------
    ! The number of the day in the year
    !---------------------------------------------
    class(t_datetime), intent(in) :: this

    call this%isLeap
    yearday = sum(daysInMonth(1:this%month - 1))
    yearday = yearday + this%day
    debug(yearday)
    call reset_DIM

  end function yearday
  !===========================================
  real(rk) function date2num(this)
    !---------------------------------------------
    ! Gives date as 'seconds from 1900-01-01 00:00:00'
    !---------------------------------------------
    integer                     :: year
    class(t_datetime), intent(in) :: this

    dbghead(date2num)
    debug(this%shortDate(.true.))

    date2num = 0
    do year = 1900, this%year - 1
      date2num = date2num + daysInYear(year) * d2s
    end do

    date2num = date2num + float(this%yearday()) * d2s + &
               float(this%hour) * h2s + float(this%minute) * min2s + &
               float(this%second) - d2s ! Why is it off by one day?
    debug(date2num)

    dbgtail(date2num)
  end function date2num
  !===========================================
  real(rk) function date_diff(start, end)
    !---------------------------------------------
    ! Gives the difference between two dates in seconds
    ! TODO: - operator?
    !---------------------------------------------
    type(t_datetime), intent(in) :: start, end
    real(rk)                     :: start_num, end_num

    start_num = start%date2num()
    end_num = end%date2num()
    date_diff = end_num - start_num

    debug(start_num)
    debug(end_num)
    debug(date_diff)

  end function date_diff
  !===========================================
  elemental real(rk) function daysInYear(year)
    !---------------------------------------------
    ! Return the number of days in year
    !---------------------------------------------
    integer, intent(in) :: year

    if (isLeapBool(year)) then
      daysInYear = 366.0d0
    else
      daysInYear = 365.0d0
    end if

  end function daysInYear
  !===========================================
  elemental integer(rk) function shortDate(this, include_time)
    !---------------------------------------------
    ! Date in short format (YYYYMMDD)
    ! If include_time is true, then include time, too (for output)
    !---------------------------------------------
    logical, intent(in)           :: include_time
    character(len=14)             :: formatout
    class(t_datetime), intent(in) :: this

    select case (include_time)
    case (.false.)
      write (formatout, '(i0.4,i0.2,i0.2)') this%year, this%month, this%day
      read (formatout, *) shortDate
    case (.true.)
      write (formatout, '(i0.4,i0.2,i0.2,i0.2,i0.2,i0.2)') this%year, this%month, this%day, this%hour, this%minute, this%second
      read (formatout, *) shortDate
    end select

  end function shortDate
  !===========================================
  elemental logical function date_gt(this, other) result(gt)
    !---------------------------------------------
    ! For comparing dates.
    ! This overloads '.gt.' or '>'
    !---------------------------------------------
    class(t_datetime), intent(in) :: this, other

    if (this%year .gt. other%year) then
      gt = .true.
    else if (this%year .lt. other%year) then
      gt = .false.
    else
      if (this%month .gt. other%month) then
        gt = .true.
      else if (this%month .lt. other%month) then
        gt = .false.
      else
        if (this%day .gt. other%day) then
          gt = .true.
        else if (this%day .lt. other%day) then
          gt = .false.
        else
          if (this%hour .gt. other%hour) then
            gt = .true.
          else if (this%hour .lt. other%hour) then
            gt = .false.
          else
            if (this%minute .gt. other%minute) then
              gt = .true.
            else if (this%minute .lt. other%minute) then
              gt = .false.
            else
              if (this%second .gt. other%second) then
                gt = .true.
              else
                gt = .false.
              end if
            end if
          end if
        end if
      end if
    end if

  end function date_gt
  !===========================================
  elemental logical function date_lt(this, other) result(lt)
    !---------------------------------------------
    ! For comparing dates.
    ! This overloads '.lt.' or '<'
    !---------------------------------------------
    class(t_datetime), intent(in) :: this, other

    lt = other > this

  end function date_lt
  !===========================================
  elemental logical function date_eq(this, other) result(eq)
    !---------------------------------------------
    ! For comparing dates
    ! This overloads '.eq.' or '=='
    !---------------------------------------------
    class(t_datetime), intent(in) :: this, other

    eq = (this%year .eq. other%year) .and. &
         (this%month .eq. other%month) .and. &
         (this%day .eq. other%day) .and. &
         (this%hour .eq. other%hour) .and. &
         (this%minute .eq. other%minute) .and. &
         (this%second .eq. other%second)

  end function date_eq
  !===========================================
  elemental logical function date_ge(this, other) result(ge)
    !---------------------------------------------
    ! You get the point.
    !---------------------------------------------
    class(t_datetime), intent(in) :: this, other

    ge = (this > other) .or. (this == other)

  end function date_ge
  !===========================================
  elemental logical function date_le(this, other) result(le)

    class(t_datetime), intent(in) :: this, other

    le = (this < other) .or. (this == other)

  end function date_le
  !===========================================
  subroutine reset_DIM()

    daysInMonth(2) = 28

    return
  end subroutine reset_DIM

end module mod_datetime
