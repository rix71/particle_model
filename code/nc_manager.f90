#define SAYLESS
#include "cppdefs.h"
module nc_manager
  !----------------------------------------------------------------
  ! Some useful netCDF subroutines
  !----------------------------------------------------------------
  use precdefs
  use netcdf
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: nc_read1d, nc_read2d, nc_read4d, nc_read_time_val, &
            nc_get_dim, nc_get_timeunit, nc_var_exists, nc_check

  !---------------------------------------------
  ! Public write variables/functions
  public :: FILLVALUE_BIG, FILLVALUE_TOPO, &
            nc_initialise, nc_add_dimension, &
            nc_add_variable, nc_add_attr, nc_write

  !---------------------------------------------
  real(rk), parameter :: FILLVALUE_TOPO = -10.0d0
  real(rk), parameter :: FILLVALUE_BIG = -9999.0d0
  !---------------------------------------------
  ! Overload the writing subroutines
  interface nc_write
    module procedure nc_write_real_1d
    module procedure nc_write_int_2d_const
    module procedure nc_write_real_2d_const
    module procedure nc_write_real_3d
    module procedure nc_write_real_4d
  end interface nc_write
  !===================================================
contains
  !===========================================
  subroutine nc_initialise(FILE_NAME)

    integer                      :: ncid
    character(len=*), intent(in) :: FILE_NAME

    FMT1, "======== Init netCDF output ========"

    call nc_check(trim(FILE_NAME), nf90_create(trim(FILE_NAME), nf90_netcdf4, ncid), "init :: create")
    call nc_check(trim(FILE_NAME), nf90_close(ncid), "init :: close")

    FMT2, trim(FILE_NAME), " initialized successfully"

    return
  end subroutine nc_initialise
  !===========================================
  subroutine nc_add_dimension(FILE_NAME, dimname, dimid, dimsize)

    integer, intent(inout)         :: dimid
    integer, intent(in), optional  :: dimsize
    integer                        :: ncid
    character(len=*), intent(in)   :: dimname
    character(len=*), intent(in)   :: FILE_NAME

    FMT1, "======== Add netCDF dimension ========"
    FMT2, "Adding dimension ", trim(dimname), " to ", trim(FILE_NAME)

    call nc_check(trim(FILE_NAME), nf90_open(trim(FILE_NAME), nf90_write, ncid), "add dim :: open")
    call nc_check(trim(FILE_NAME), nf90_redef(ncid), "add dim :: redef mode")
    if (present(dimsize)) then
      call nc_check(trim(FILE_NAME), nf90_def_dim(ncid, trim(dimname), dimsize, dimid), "add dim :: def "//trim(dimname))
    else
      call nc_check(trim(FILE_NAME), nf90_def_dim(ncid, trim(dimname), nf90_unlimited, dimid), "add dim :: def "//trim(dimname))
    end if
    call nc_check(trim(FILE_NAME), nf90_enddef(ncid), "add dim :: end def")
    call nc_check(trim(FILE_NAME), nf90_close(ncid), "add dim :: close")

    FMT2, trim(dimname), " added successfully"

    return
  end subroutine nc_add_dimension
  !===========================================
  subroutine nc_add_variable(FILE_NAME, varname, dType, nDims, dimids, missing_val)
    !---------------------------------------------
    ! Add variables to output file
    !---------------------------------------------

    character(len=*), intent(in)   :: FILE_NAME
    character(len=*), intent(in)   :: varname
    character(len=*), intent(in)   :: dType
    integer, intent(in)            :: nDims
    integer, intent(in)            :: dimids(nDims)
    integer                        :: ncid, varid
    real(rk), intent(in), optional :: missing_val

    FMT1, "======== Add netCDF variable ========"
    FMT2, "Adding variable ", trim(varname), " to ", trim(FILE_NAME)

    call nc_check(trim(FILE_NAME), nf90_open(trim(FILE_NAME), nf90_write, ncid), "add var :: open")
    call nc_check(trim(FILE_NAME), nf90_redef(ncid), "add var :: redef mode")
    select case (dType)
    case ('float')
      call nc_check(trim(FILE_NAME), nf90_def_var(ncid, varname, nf90_double, dimids, varid), "add var :: def "//trim(varname))
      if (present(missing_val)) then
        call nc_check(trim(FILE_NAME), nf90_put_att(ncid, varid, "missing_value", missing_val), "add var :: def "//trim(varname))
      end if
    case ('int')
      call nc_check(trim(FILE_NAME), nf90_def_var(ncid, varname, nf90_int, dimids, varid), "add var :: def "//trim(varname))
      if (present(missing_val)) then
      call nc_check(trim(FILE_NAME), nf90_put_att(ncid, varid, "missing_value", int(missing_val)), "add var :: def "//trim(varname))
      end if
    end select
    call nc_check(trim(FILE_NAME), nf90_enddef(ncid), "add var :: end def")
    call nc_check(trim(FILE_NAME), nf90_close(ncid), "add var :: close")

    FMT2, trim(varname), " added successfully"

    return
  end subroutine nc_add_variable
  !===========================================
  subroutine nc_add_attr(FILE_NAME, varname, attrname, attrval)
    !---------------------------------------------
    ! Add attribute to variable
    !---------------------------------------------

    character(len=*), intent(in)   :: FILE_NAME
    character(len=*), intent(in)   :: varname, attrname, attrval
    integer                        :: ncid, varid

    FMT1, "======== Add netCDF attribute ========"
    FMT2, "Adding attribute ", trim(attrname), " to ", trim(FILE_NAME)

    call nc_check(trim(FILE_NAME), nf90_open(trim(FILE_NAME), nf90_write, ncid), "add attr :: open")
    call nc_check(trim(FILE_NAME), nf90_redef(ncid), "add attr :: redef mode")
    call nc_check(trim(FILE_NAME), nf90_inq_varid(ncid, trim(varname), varid), "add attr :: inq varid")
    call nc_check(trim(FILE_NAME), nf90_put_att(ncid, varid, trim(attrname), trim(attrval)), "add attr :: put "//trim(attrname))
    call nc_check(trim(FILE_NAME), nf90_enddef(ncid), "add attr :: end def")
    call nc_check(trim(FILE_NAME), nf90_close(ncid), "add attr :: close")

    FMT2, trim(attrname), " added successfully"

    return
  end subroutine nc_add_attr
  !===========================================
  subroutine nc_write_real_1d(FILE_NAME, datain, varname, nvals)
    !---------------------------------------------
    ! Write 1D real data with no time axis (output will be 1D)
    !---------------------------------------------

    integer, parameter           :: nDims = 1
    integer                      :: ncid, varid, start(nDims), count(nDims)
    integer, intent(in)          :: nvals
    character(len=*), intent(in) :: FILE_NAME
    character(len=*), intent(in) :: varname
    real(rk), intent(in)         :: datain(nvals)

    FMT1, "======== Write netCDF variable ========"
    FMT2, "Writing variable ", trim(varname), " to ", trim(FILE_NAME)

    call nc_check(trim(FILE_NAME), nf90_open(trim(FILE_NAME), nf90_write, ncid), "write :: open")
    call nc_check(trim(FILE_NAME), nf90_inq_varid(ncid, trim(varname), varid), "write :: inq varid")
    start = (/1/)
    count = (/nvals/)
    call nc_check(trim(FILE_NAME), nf90_put_var(ncid, varid, datain, start=start, count=count), "write :: put var")
    call nc_check(trim(FILE_NAME), nf90_close(ncid), "write :: close")

    FMT2, trim(varname), " written successfully"

    return
  end subroutine nc_write_real_1d
  !===========================================
  subroutine nc_write_int_2d_const(FILE_NAME, datain, varname, nx, ny)
    !---------------------------------------------
    ! Write 2D integer data with no time axis (output will be 2D)
    !---------------------------------------------

    integer, intent(in)          :: nx, ny
    integer, intent(in)          :: datain(nx, ny)
    integer, parameter           :: nDims = 2
    integer                      :: ncid, varid, start(nDims), count(nDims)
    character(len=*), intent(in) :: FILE_NAME
    character(len=*), intent(in) :: varname

    FMT1, "======== Write netCDF variable ========"
    FMT2, "Writing variable ", trim(varname), " to ", trim(FILE_NAME)

    call nc_check(trim(FILE_NAME), nf90_open(trim(FILE_NAME), nf90_write, ncid), "write :: open")
    call nc_check(trim(FILE_NAME), nf90_inq_varid(ncid, trim(varname), varid), "write :: inq varid")
    start = (/1, 1/)
    count = (/nx, ny/)
    call nc_check(trim(FILE_NAME), nf90_put_var(ncid, varid, datain, start=start, count=count), "write :: put var")
    call nc_check(trim(FILE_NAME), nf90_close(ncid), "write :: close")

    FMT2, trim(varname), " written successfully"

    return
  end subroutine nc_write_int_2d_const
  !===========================================
  subroutine nc_write_real_2d_const(FILE_NAME, datain, varname, nx, ny)
    !---------------------------------------------
    ! Write 2D real data with no time axis (output will be 2D)
    !---------------------------------------------

    integer, parameter           :: nDims = 2
    integer                      :: ncid, varid, start(nDims), count(nDims)
    integer, intent(in)          :: nx, ny
    character(len=*), intent(in) :: FILE_NAME
    character(len=*), intent(in) :: varname
    real(rk), intent(in)         :: datain(nx, ny)

    FMT1, "======== Write netCDF variable ========"
    FMT2, "Writing variable ", trim(varname), " to ", trim(FILE_NAME)

    call nc_check(trim(FILE_NAME), nf90_open(trim(FILE_NAME), nf90_write, ncid), "write :: open")
    call nc_check(trim(FILE_NAME), nf90_inq_varid(ncid, trim(varname), varid), "write :: inq varid")
    start = (/1, 1/)
    count = (/nx, ny/)
    call nc_check(trim(FILE_NAME), nf90_put_var(ncid, varid, datain, start=start, count=count), "write :: put var")
    call nc_check(trim(FILE_NAME), nf90_close(ncid), "write :: close")

    FMT2, trim(varname), " written successfully"

    return
  end subroutine nc_write_real_2d_const
  !===========================================
  subroutine nc_write_real_3d(FILE_NAME, datain, varname, nx, ny, itime)
    !---------------------------------------------
    ! Write 2D real data to itime time step (output will be 3D)
    !---------------------------------------------

    integer, intent(in)          :: nx, ny, itime
    integer, parameter           :: nDims = 3
    integer                      :: ncid, varid, start(nDims), count(nDims)
    character(len=*), intent(in) :: FILE_NAME
    character(len=*), intent(in) :: varname
    real(rk), intent(in)         :: datain(nx, ny)

    FMT1, "======== Write netCDF variable ========"
    FMT2, "Writing variable ", trim(varname), " to ", trim(FILE_NAME)

    call nc_check(trim(FILE_NAME), nf90_open(trim(FILE_NAME), nf90_write, ncid), "write :: open")
    call nc_check(trim(FILE_NAME), nf90_inq_varid(ncid, trim(varname), varid), "write :: inq varid")
    start = [1, 1, itime]
    count = [nx, ny, 1]
    call nc_check(trim(FILE_NAME), nf90_put_var(ncid, varid, datain, start=start, count=count), "write :: put var")
    call nc_check(trim(FILE_NAME), nf90_close(ncid), "write :: close")

    FMT2, trim(varname), " written successfully"

    return
  end subroutine nc_write_real_3d
  !===========================================
  subroutine nc_write_real_4d(FILE_NAME, datain, varname, nx, ny, nz, itime)
    !---------------------------------------------
    ! Write 3D real data to itime time step (output will be 4D)
    !---------------------------------------------

    integer, intent(in)          :: nx, ny, nz, itime
    integer, parameter           :: nDims = 4
    integer                      :: ncid, varid, start(nDims), count(nDims)
    character(len=*), intent(in) :: FILE_NAME
    character(len=*), intent(in) :: varname
    real(rk), intent(in)         :: datain(nx, ny, nz)

    FMT1, "======== Write netCDF variable ========"
    FMT2, "Writing variable ", trim(varname), " to ", trim(FILE_NAME)

    call nc_check(trim(FILE_NAME), nf90_open(trim(FILE_NAME), nf90_write, ncid), "write :: open")
    call nc_check(trim(FILE_NAME), nf90_inq_varid(ncid, trim(varname), varid), "write :: inq varid")
    start = (/1, 1, 1, itime/)
    count = (/nx, ny, nz, 1/)
    call nc_check(trim(FILE_NAME), nf90_put_var(ncid, varid, datain, start=start, count=count), "write :: put var")
    call nc_check(trim(FILE_NAME), nf90_close(ncid), "write :: close")

    FMT2, trim(varname), " written successfully"

    return
  end subroutine nc_write_real_4d
  !===========================================
  subroutine nc_read1d(fname, vname, nvals, dataout)

    integer               :: nvals
    integer               :: ncid, varid
    character(len=*)      :: fname
    character(len=*)      :: vname
    real(rk), intent(out) :: dataout(nvals)

    call nc_check(trim(fname), nf90_open(fname, nf90_nowrite, ncid), "nc_read1d :: open")
    call nc_check(trim(fname), nf90_inq_varid(ncid, vname, varid), "nc_read1d :: inq_varid "//trim(vname))
    call nc_check(trim(fname), nf90_get_var(ncid, varid, dataout), "nc_read1d :: get_var "//trim(vname))
    call nc_check(trim(fname), nf90_close(ncid), "nc_read1d :: close")

    return
  end subroutine
  !===========================================
  subroutine nc_read2d(fname, vname, nx, ny, dataout)

    integer              :: nx, ny
    integer              :: ncid, varid
    character(len=*)     :: fname
    character(len=*)     :: vname
    real(rk), intent(out) :: dataout(nx, ny)

    call nc_check(trim(fname), nf90_open(fname, nf90_nowrite, ncid), "nc_read2d :: open")
    call nc_check(trim(fname), nf90_inq_varid(ncid, vname, varid), "nc_read2d :: inq_varid "//trim(vname))
    call nc_check(trim(fname), nf90_get_var(ncid, varid, dataout), "nc_read2d :: get_var "//trim(vname))
    call nc_check(trim(fname), nf90_close(ncid), "nc_read2d :: close")

    return
  end subroutine nc_read2d
  !===========================================
  subroutine nc_read4d(fname, vname, start, count, dataout)

    integer               :: ncid, varid
    integer, dimension(4) :: start, count
    character(len=*)      :: fname
    character(len=*)      :: vname
    real(rk), intent(out) :: dataout(count(1), count(2), count(3), count(4))

    call nc_check(trim(fname), nf90_open(fname, nf90_nowrite, ncid), "nc_read4d :: open")
    call nc_check(trim(fname), nf90_inq_varid(ncid, vname, varid), "nc_read4d :: inq_varid "//trim(vname))
    call nc_check(trim(fname), nf90_get_var(ncid, varid, dataout, start=start, &
                                            count=count), "nc_read4d :: get_var "//trim(vname))
    call nc_check(trim(fname), nf90_close(ncid), "nc_read4d :: close")

    return
  end subroutine
  !===========================================
  subroutine nc_read_time_val(fname, n, timeval)

    integer                      :: ncid, varid
    integer, intent(in)          :: n
    character(len=*), intent(in) :: fname
    real(rk), intent(out)        :: timeval
    real(rk)                     :: tmpval(1)

    call nc_check(trim(fname), nf90_open(fname, nf90_nowrite, ncid), "nc_read_time_val :: open")
    call nc_check(trim(fname), nf90_inq_varid(ncid, 'time', varid), "nc_read_time_val :: inq_var_id 'time'")
    call nc_check(trim(fname), nf90_get_var(ncid, varid, tmpval, start=[n], count=[1]), "nc_read_time_val :: get_var 'time'")
    call nc_check(trim(fname), nf90_close(ncid), "nc_read_time_val :: close")

    timeval = tmpval(1)

    return
  end subroutine nc_read_time_val
  !===========================================
  subroutine nc_get_dim(fname, dname, ndim)

    integer, intent(out)         :: ndim
    integer                      :: ncid, dimid
    character(len=*), intent(in) :: fname
    character(len=*), intent(in) :: dname

    call nc_check(trim(fname), nf90_open(fname, nf90_nowrite, ncid), "get_dim :: open")
    call nc_check(trim(fname), nf90_inq_dimid(ncid, dname, dimid), 'get_dim :: inq_dim_id '//trim(dname))
    call nc_check(trim(fname), nf90_inquire_dimension(ncid, dimid, len=ndim), 'get_dim :: inq_dim '//trim(dname))
    call nc_check(trim(fname), nf90_close(ncid), 'get_dim :: close')

    return
  end subroutine nc_get_dim
  !===========================================
  subroutine nc_get_timeunit(fname, timeunit)

    integer                       :: ncid, varid
    character(len=*), intent(in)  :: fname
    character(len=*), intent(out) :: timeunit

    call nc_check(trim(fname), nf90_open(fname, nf90_nowrite, ncid), "get_timeunit :: open")
    call nc_check(trim(fname), nf90_inq_varid(ncid, "time", varid), "get_timeunit :: inq_var_id 'time'")
    call nc_check(trim(fname), nf90_get_att(ncid, varid, 'units', timeunit), "get_timeunit :: get_attr 'units'")
    call nc_check(trim(fname), nf90_close(ncid), 'get_timeunit :: close')

    return
  end subroutine nc_get_timeunit
  !===========================================
  logical function nc_var_exists(fname, vname)

    character(len=*), intent(in) :: fname
    character(len=*), intent(in) :: vname
    integer                      :: ncid, varid

    nc_var_exists = .true.
    call nc_check(trim(fname), nf90_open(fname, nf90_nowrite, ncid), "nc_var_exists :: open")
    if (nf90_inq_varid(ncid, trim(vname), varid) == nf90_enotvar) nc_var_exists = .false.
    call nc_check(trim(fname), nf90_close(ncid), 'nc_var_exists :: close')

  end function nc_var_exists
  !===========================================
  subroutine nc_check(fname, status, code)

    integer, intent(in) :: status
    character(len=*)    :: fname, code

    if (status /= nf90_noerr) then
      ERROR, "NETCDF: Stopped at "//trim(code)//" with ", status
      ERROR, trim(nf90_strerror(status))
      ERROR, "NETCDF: File name: "//trim(fname)
      stop
    end if

  end subroutine nc_check

end module nc_manager
