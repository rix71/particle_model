#include "cppdefs.h"
module mod_fieldset
  use mod_precdefs
  use mod_errors
  use mod_datetime
  use nc_manager
  use mod_field, only: t_field
  use mod_domain, only: t_domain
  use mod_list, only: t_list
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_fieldset
  !---------------------------------------------
  character(len=*), parameter     :: dirinfile = 'dirlist.dat' ! List the dirs into this file
  !---------------------------------------------
  type t_fieldset
    private
    !---------------------------------------------
    type(t_list)                    :: fields
    integer, public                 :: num_fields = 0
    integer, public                 :: nx, ny, nz
    real(rk)                        :: nc_timestep
    ! real(rk), allocatable         :: bathymetry(:, :)
    type(t_domain), pointer, public :: domain
    type(t_datetime)                :: date_t1, date_t2
    !---------------------------------------------
    ! Velocity component variables
    integer :: u_idx = 0
    integer :: v_idx = 0
    integer, allocatable :: u_mask(:, :)
    integer, allocatable :: v_mask(:, :)
    !---------------------------------------------
    ! Z axis variables
    integer      :: zax_style = DEPTH_VALUES
    integer      :: zax_idx = 0
    !---------------------------------------------
    ! netCDF read variables
    logical                   :: read_first = .true.
    integer                   :: read_idx
    integer                   :: read_idx_increment = 1
    character(len=LEN_CHAR_L) :: current_path   ! Path to current data
    integer                   :: current_ntimes ! Current timestep index
    integer                   :: dirlist_idx
    type(t_datetime)          :: next_read_dt
    !---------------------------------------------
    ! Paths/directories
    character(len=LEN_CHAR_L)              :: PATH, PMAPFILE           ! Path to field data and processor map
    character(len=LEN_CHAR_L)              :: file_prefix, file_suffix ! What comes before and after the proc. number?
    character(len=LEN_CHAR_L), allocatable :: filelist(:)              ! List of files
    integer, allocatable                   :: dirlist(:)               ! List of directories (one or the other will be allocated)
    integer                                :: nentries                 ! Number of directories or files
    !---------------------------------------------
    ! Subdomain variables
    logical              :: has_subdomains = .false.
    integer              :: nproc       ! Number of subdomains
    integer              :: nxp, nyp    ! Size of subdomain
    integer, allocatable :: pmap(:, :)  ! Offsets of subdomains
    integer, allocatable :: pmask(:, :) ! Proc. mask

  contains
    private
    procedure, public :: add_field
    procedure, public :: list_fields
    procedure, public :: get_field_index
    procedure, public :: set_u_component
    procedure, public :: set_v_component
    procedure         :: get_u_comp_mask
    procedure         :: get_v_comp_mask
    procedure, public :: set_zax
    procedure, public :: get_zax
    procedure, public :: get_gradient
    procedure, public :: set_start_time
    procedure, public :: set_simulation_timestep
    procedure, public :: get_folder, get_file
    procedure, public :: get_time
    procedure, public :: find_folder, find_file
    procedure, public :: get_pmap
    procedure         :: get_indices_vertical
    procedure, public :: search_indices
    procedure, public :: sealevel
    procedure         :: init_proc_mask
    procedure         :: init_dirlist
    generic, public   :: set => set_value_key, set_value_idx
    procedure         :: set_value_key, set_value_idx
    generic, public   :: get => get_value_key_real_idx, get_value_key_int_idx, get_value_idx_real_idx, get_value_idx_int_idx
    procedure         :: get_value_key_real_idx, get_value_key_int_idx, get_value_idx_real_idx, get_value_idx_int_idx
    procedure         :: read_field
    procedure         :: read_field_subdomains
    procedure, public :: update
    procedure, public :: read_first_timesteps

  end type t_fieldset
  !---------------------------------------------
  interface t_fieldset
    module procedure :: ctor_fieldset
  end interface t_fieldset
  !---------------------------------------------
  integer :: ierr
  !===================================================
contains
  !===========================================
  type(t_fieldset) function ctor_fieldset(nx, ny, nz, file_prefix, file_suffix, domain, path, pmap) result(f)
#undef PROC0
#define PROC0 "/"//trim(f%file_prefix)//"0000"//trim(f%file_suffix)//".nc"
    integer, intent(in) :: nx, ny, nz
    character(*), intent(in) :: file_prefix, file_suffix
    type(t_domain), target, intent(in) :: domain
    character(*), optional, intent(in) :: path, pmap
    character(len=LEN_CHAR_L) :: init_path
    real(rk) :: t1, t2

    f%nx = nx
    f%ny = ny
    f%nz = nz
    ! allocate (f%domain(nx, ny))
    f%domain => domain
    if (present(path)) f%PATH = path
    if (present(pmap)) then
      f%PMAPFILE = pmap
      f%has_subdomains = .true.
      call f%init_proc_mask()
    end if
    f%file_prefix = file_prefix
    f%file_suffix = file_suffix
    call f%init_dirlist()

    select case (f%has_subdomains)
    case (.true.)
      init_path = f%get_folder(1)
      t1 = nc_read_time_val(trim(init_path)//PROC0, 1)
      t2 = nc_read_time_val(trim(init_path)//PROC0, 2)
    case (.false.)
      init_path = f%get_file(1)
      t1 = nc_read_time_val(trim(init_path), 1)
      t2 = nc_read_time_val(trim(init_path), 2)
    end select
    f%nc_timestep = t2 - t1
    if (f%nc_timestep <= ZERO) then

      ERROR, "Dataset time step <= 0: ", f%nc_timestep, " [s]"
      call throw_error("fieldset :: fieldset", "Bad time step")
    end if

    return
#undef PROC0
#define PROC0 "/"//trim(this%file_prefix)//"0000"//trim(this%file_suffix)//".nc"
  end function ctor_fieldset
  !===========================================
  subroutine add_field(this, field_name, nc_varname, is_2d)
    class(t_fieldset), intent(inout) :: this
    character(*), intent(in)         :: field_name, nc_varname
    logical, optional, intent(in)    :: is_2d
    logical :: fld_2d
#ifdef DEBUG
    integer, save :: NFIELD = 0

    NFIELD = NFIELD + 1

#endif

    fld_2d = .false.
    if (present(is_2d)) then
      fld_2d = is_2d
    end if

    select case (fld_2d)
    case (.false.)

      call this%fields%add_node(field_name, &
                                t_field(this%nx, this%ny, this%nz, timestep=this%nc_timestep, nc_varname=nc_varname))
    case (.true.)

      call this%fields%add_node(field_name, &
                                t_field(this%nx, this%ny, timestep=this%nc_timestep, nc_varname=nc_varname))
    end select

    this%num_fields = this%num_fields + 1

    return
  end subroutine add_field
  !===========================================
  subroutine set_value_key(this, field_name, data)
    class(t_fieldset), intent(inout) :: this
    character(*), intent(in)         :: field_name
    real(rk), intent(in)             :: data(this%nx, this%ny, this%nz)
    type(t_field), pointer           :: p_field

    call this%fields%get_item(field_name, p_field)
    call p_field%set(data)

  end subroutine set_value_key
  !===========================================
  subroutine set_value_idx(this, idx, data)
    class(t_fieldset), intent(inout) :: this
    integer, intent(in)              :: idx
    real(rk), intent(in)             :: data(this%nx, this%ny, this%nz)
    type(t_field), pointer           :: p_field

    call this%fields%get_item(idx, p_field)
    call p_field%set(data)

  end subroutine set_value_idx
  !===========================================
  real(rk) function get_value_key_real_idx(this, field_name, t, i, j, k) result(res)
    class(t_fieldset), intent(in)  :: this
    character(*), intent(in)       :: field_name
    real(rk), intent(in)           :: t, i, j
    real(rk), optional, intent(in) :: k
    type(t_field), pointer         :: p_field

    call this%fields%get_item(field_name, p_field)
    if (present(k)) then

      call p_field%get(t=t, x=i, y=j, z=k, &
                       res=res)
    else

      call p_field%get(t=t, x=i, y=j, &
                       res=res)
    end if

    return
  end function get_value_key_real_idx
  !===========================================
  real(rk) function get_value_key_int_idx(this, field_name, t, i, j, k) result(res)
    class(t_fieldset), intent(in)  :: this
    character(*), intent(in)       :: field_name
    real(rk), intent(in)           :: t
    integer, intent(in)            :: i, j
    integer, optional, intent(in)  :: k
    type(t_field), pointer         :: p_field

    call this%fields%get_item(field_name, p_field)
    if (present(k)) then

      call p_field%get(t=t, i=i, j=j, k=k, res=res)
    else

      call p_field%get(t=t, i=i, j=j, res=res)
    end if

    return
  end function get_value_key_int_idx
  !===========================================
  real(rk) function get_value_idx_real_idx(this, idx, t, i, j, k) result(res)
    class(t_fieldset), intent(in)  :: this
    integer, intent(in)            :: idx
    real(rk), intent(in)           :: t, i, j
    real(rk), optional, intent(in) :: k
    type(t_field), pointer         :: p_field

    call this%fields%get_item(idx, p_field)
    if (present(k)) then
      call p_field%get(t, i, j, z=k, &
                       res=res)
    else
      call p_field%get(t, i, j, &
                       res=res)
    end if

  end function get_value_idx_real_idx
  !===========================================
  real(rk) function get_value_idx_int_idx(this, idx, t, i, j, k) result(res)
    class(t_fieldset), intent(in)  :: this
    integer, intent(in)            :: idx
    real(rk), intent(in)           :: t
    integer, intent(in)            :: i, j
    integer, optional, intent(in)  :: k
    type(t_field), pointer         :: p_field

    call this%fields%get_item(idx, p_field)
    if (present(k)) then
      call p_field%get(t, i=i, j=j, k=k, res=res)
    else
      call p_field%get(t, i=i, j=j, res=res)
    end if

  end function get_value_idx_int_idx
  !===========================================
  subroutine list_fields(this)
    class(t_fieldset), intent(in) :: this

    call this%fields%get_info()

  end subroutine list_fields
  !===========================================
  integer function get_field_index(this, field_name) result(res)
    class(t_fieldset), intent(in) :: this
    character(*), intent(in)      :: field_name

    res = this%fields%node_loc(field_name)

    return
  end function get_field_index
  !===========================================
  subroutine get_u_comp_mask(this)
    !---------------------------------------------
    ! Make the velocity u component (similarly for v comp. in get_v_comp_mask)
    ! to deal with boundaries.
    ! SEA will be 0, BEACH areas will be 1 or -1 depending on where the land is.
    ! Adjacent land cells will have the same velocity tangential to land (0 gradient).
    ! Normal (cross-boundary) components will be 0.
    !---------------------------------------------
    class(t_fieldset), intent(inout) :: this
    integer                          :: i, j

    allocate (this%u_mask(this%nx, this%ny), stat=ierr)
    if (ierr .ne. 0) call throw_error("fieldset :: get_u_comp_mask", "Could not allocate u component mask")
    this%u_mask = 0
    do i = 1, this%nx
      do j = 1, this%ny
        if (this%domain%get_seamask(i, j) == BEACH) then
          if (this%domain%get_seamask(i, j + 1) == LAND) then
            this%u_mask(i, j) = -1
          else if (this%domain%get_seamask(i, j - 1) == LAND) then
            this%u_mask(i, j) = 1
          end if
        end if
      end do
    end do

  end subroutine get_u_comp_mask
  !===========================================
  subroutine get_v_comp_mask(this)
    class(t_fieldset), intent(inout) :: this
    integer                          :: i, j

    allocate (this%v_mask(this%nx, this%ny), stat=ierr)
    if (ierr .ne. 0) call throw_error("fieldset :: get_v_comp_mask", "Could not allocate v component mask")
    this%v_mask = 0
    do i = 1, this%nx
      do j = 1, this%ny
        if (this%domain%get_seamask(i, j) == BEACH) then
          if (this%domain%get_seamask(i + 1, j) == LAND) then
            this%v_mask(i, j) = -1
          else if (this%domain%get_seamask(i - 1, j) == LAND) then
            this%v_mask(i, j) = 1
          end if
        end if
      end do
    end do

  end subroutine get_v_comp_mask
  !===========================================
  subroutine set_u_component(this, u_comp_name)
    !---------------------------------------------
    ! Set the index of the velocity u component and create
    ! the mask (similarly for v component in set_v_component).
    !---------------------------------------------

    class(t_fieldset), intent(inout) :: this
    character(*), intent(in)         :: u_comp_name

    this%u_idx = this%fields%node_loc(trim(u_comp_name))
    if (this%u_idx < 1) call throw_error("fieldset :: set_u_component", "Did not find "//trim(u_comp_name)//" in fieldset")
    call this%get_u_comp_mask()

  end subroutine set_u_component
  !===========================================
  subroutine set_v_component(this, v_comp_name)
    class(t_fieldset), intent(inout) :: this
    character(*), intent(in)         :: v_comp_name

    this%v_idx = this%fields%node_loc(trim(v_comp_name))
    if (this%v_idx < 1) call throw_error("fieldset :: set_v_component", "Did not find "//trim(v_comp_name)//" in fieldset")
    call this%get_v_comp_mask()

  end subroutine set_v_component
  !===========================================
  subroutine set_zax(this, zax_name, zax_style)
    class(t_fieldset), intent(inout) :: this
    character(*), intent(in)         :: zax_name
    integer, intent(in)              :: zax_style

    this%zax_style = zax_style
    this%zax_idx = this%fields%node_loc(trim(zax_name))

    if (this%zax_idx < 1) call throw_error("fieldset :: set_zax", "Did not find "//trim(zax_name)//" in fieldset")

    return
  end subroutine set_zax
  !===========================================
  subroutine set_start_time(this, date)
    class(t_fieldset), intent(inout) :: this
    type(t_datetime), intent(in)     :: date

    if (this%has_subdomains) then
      call this%find_folder(date, this%current_path, this%dirlist_idx)
      call nc_get_dim(trim(this%current_path)//PROC0, 'time', this%current_ntimes)
      if (nc_read_time_val(trim(this%current_path)//PROC0, 1) == ZERO) then
        ! If RefTime is the same as the first nc_timestep, then 1 must be added
        this%read_idx = int(date_diff(datetime_from_netcdf(trim(this%current_path)//PROC0), date) / this%nc_timestep) + 1
      else
        this%read_idx = int(date_diff(datetime_from_netcdf(trim(this%current_path)//PROC0), date) / this%nc_timestep)
      end if
    else
      call this%find_file(date, this%current_path, this%dirlist_idx)
      call nc_get_dim(trim(this%current_path), 'time', this%current_ntimes)
      if (nc_read_time_val(trim(this%current_path), 1) == ZERO) then
        this%read_idx = int(date_diff(datetime_from_netcdf(trim(this%current_path)), date) / this%nc_timestep) + 1
      else
        this%read_idx = int(date_diff(datetime_from_netcdf(trim(this%current_path)), date) / this%nc_timestep)
      end if
    end if
    this%next_read_dt = date
    this%date_t1 = date
    this%date_t2 = date%nextDate(this%nc_timestep)

    return
  end subroutine set_start_time
  !===========================================
  subroutine set_simulation_timestep(this, dt)
    class(t_fieldset), intent(inout) :: this
    real(rk), intent(in)             :: dt

    if (dt <= this%nc_timestep) then
      this%read_idx_increment = 1
    else
      if (mod(dt, this%nc_timestep) .ne. 0) then
        ERROR, "Time step ", dt, " [s] not divisible by netCDF time step ", this%nc_timestep, " [s]"
        call throw_error("fieldset :: set_simulation_timestep", "Time step should be divisible by netCDF time step")
      end if
      this%read_idx_increment = int(dt / this%nc_timestep)
    end if

    return
  end subroutine set_simulation_timestep
  !===========================================
  real(rk) function get_time(this, date) result(res)
    class(t_fieldset), intent(in) :: this
    type(t_datetime), intent(in)  :: date

    res = date_diff(this%date_t1, date)

  end function get_time
  !===========================================
  character(len=LEN_CHAR_L) function get_folder(this, idx) result(res)
    class(t_fieldset), intent(in) :: this
    integer, intent(in)           :: idx

    if (idx < 1) call throw_error("fieldset :: get_folder", "Index out of bounds (idx < 1)")
    if (idx > this%nentries) call throw_error("fieldset :: get_folder", "Index out of bounds (idx > number of folders)")
    write (res, '(a,i0.8)') trim(this%PATH)//'/', this%dirlist(idx)

    return
  end function get_folder
  !===========================================
  character(len=LEN_CHAR_L) function get_file(this, idx) result(res)
    class(t_fieldset), intent(in) :: this
    integer, intent(in)           :: idx

    if (idx < 1) call throw_error("fieldset :: get_file", "Index out of bounds (idx < 1)")
    if (idx > this%nentries) call throw_error("fieldset :: get_file", "Index out of bounds (idx > number of files)")
    write (res, '(a)') trim(this%PATH)//'/'//trim(this%filelist(idx))

    return
  end function get_file
  !===========================================
  subroutine init_dirlist(this)
    !---------------------------------------------
    ! Get a list of directories or files that contain the data
    ! The files should be named so that ls command would give them
    ! in the right order. Including the date in the file name should be enough.
    ! TODO: Alternatively, if all the files are in one directory,
    !       it should list all the files in this%PATH.
    !       (separate routine e.g., init_filelist ?)
    ! TODO (later, probably never): sort the files somehow so all of this
    !       would not depend on ls getting it right.
    !---------------------------------------------
    class(t_fieldset), intent(inout) :: this
    logical                          :: dirlist_exists
    integer                          :: idir

    FMT1, "======== Init dirlist ========"
    inquire (file=trim(dirinfile), exist=dirlist_exists)
    if (.not. dirlist_exists) then
#ifndef NOSYSCALLS
      select case (this%has_subdomains)
      case (.true.)
        FMT2, "Getting directory list..."
        call system('( ls -d '//trim(this%PATH)//'/*/ | wc -l && ls -d ' &
                    //trim(this%PATH)//'/*/ | xargs -n 1 basename ) > ' &
                    //trim(dirinfile), status=ierr)
        if (ierr .ne. 0) call throw_error("fieldset :: init_dirlist", "Could not get directory list from "//trim(this%PATH), ierr)
      case (.false.)
        FMT2, "Getting file list..."
        call system('( ls '//trim(this%PATH)//'/'//trim(this%file_prefix)//'*' &
                    //trim(this%file_suffix)//'*.nc | wc -l && ls '//trim(this%PATH) &
                    //'/'//trim(this%file_prefix)//'*'//trim(this%file_suffix) &
                    //'*.nc | xargs -n 1 basename ) > '//trim(dirinfile), status=ierr)
        if (ierr .ne. 0) call throw_error("fieldset :: init_dirlist", "Could not get file list from "//trim(this%PATH), ierr)
      end select
#else
      call throw_error("fieldset :: init_dirlist", "File/directory list ("//trim(dirinfile)//") does not exist!")
#endif
    end if

    open (DIRFILE, file=trim(dirinfile), action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error("fieldset :: init_dirlist", "Failed to open "//trim(dirinfile), ierr)
    read (DIRFILE, *) this%nentries
    FMT2, "Found ", this%nentries, "files/directories in "//trim(this%PATH)
    select case (this%has_subdomains)
    case (.true.)
      allocate (this%dirlist(this%nentries))
      do idir = 1, this%nentries
        read (DIRFILE, *) this%dirlist(idir)
      end do
    case (.false.)
      allocate (this%filelist(this%nentries))
      do idir = 1, this%nentries
        read (DIRFILE, '(a)') this%filelist(idir)
      end do
    end select
    close (DIRFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error("fieldset :: init_dirlist", "Failed to close "//trim(dirinfile), ierr)

    FMT2, "Dirlist initialized"

  end subroutine init_dirlist
  !===========================================
  subroutine init_proc_mask(this)
    !---------------------------------------------
    ! This maps the pieces of GETM data using par_setup.
    ! TODO: Since pmask takes up a lot of unnecessary space,
    !       maybe this should only be called when par_setup is wanted
    !      (either a compilation flag or "if (parallel)" or something in init_model)
    ! EDIT: This is called from init_model only if this%has_subdomains=.true.
    !---------------------------------------------
    class(t_fieldset), intent(inout) :: this
    integer                          :: imax, jmax ! Total size of the domain
    integer                          :: ioff, joff ! Subdomain offset
    integer                          :: pnum, iend, jend
    integer                          :: i, j, iproc

    FMT1, "======== Init subdomains ========"

    open (PROCFILE, file=trim(this%PMAPFILE), action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error("fieldset :: init_proc_mask", "Failed to open "//trim(this%PMAPFILE), ierr)
    read (PROCFILE, *) this%nproc
    allocate (this%pmap(this%nproc, 2))
    read (PROCFILE, *) this%nxp, this%nyp, imax, jmax
    if ((imax .ne. this%nx) .or. (jmax .ne. this%ny)) then
      ERROR, "(nx = ", this%nx, "imax = ", imax, "ny = ", this%ny, "jmax = ", jmax, ")"
      call throw_error("fieldset :: init_proc_mask", "Given domain size does not match one in parallel setup file!")
    else
      FMT2, "Number of subdomains: ", this%nproc, " (", this%nxp, " x ", this%nyp, ")"
    end if
    allocate (this%pmask(imax, jmax))
    do i = 1, imax
      do j = 1, jmax
        this%pmask(i, j) = -10
      end do
    end do
    iproc = 0
    do while (iproc .lt. this%nproc)
      read (PROCFILE, *) pnum, ioff, joff
      i = ioff + 1; j = joff + 1
      iproc = iproc + 1
      this%pmap(iproc, 1) = ioff; this%pmap(iproc, 2) = joff
      iend = this%nxp; jend = this%nyp
      if (i .lt. 1) then
        !---------------------------------------------
        ! If i is less than 1, then ioff must be negative.
        ! Therefore this processor actually covers less area,
        ! i.e. iend must be reduced ioff amount
        iend = this%nxp + ioff ! Same as nxp - abs(ioff)
        i = 1
      end if
      if (j .lt. 1) then
        !---------------------------------------------
        ! Same logic here
        jend = this%nyp + joff
        j = 1
      end if
      this%pmask(i:i + iend, j:j + jend) = pnum
    end do
    close (PROCFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error("fieldset :: init_proc_mask", "Failed to close "//trim(this%PMAPFILE), ierr)

    FMT2, "Processor mask initialised"

    return
  end subroutine init_proc_mask
  !===========================================
  subroutine get_pmap(this, pmapout)
    !---------------------------------------------
    ! Just in case I want to use pmap somewhere else
    ! (it's private otherwise)
    !---------------------------------------------
    class(t_fieldset), intent(in) :: this
    integer, intent(inout)        :: pmapout(this%nproc, 2)

    pmapout = this%pmap

    return
  end subroutine get_pmap
  !===========================================
  subroutine find_folder(this, date, thedir, folder_idx)
    !---------------------------------------------
    ! TODO: might have to check if nc files start at time
    ! [date] 00:00:00 or [date] 00:00:10 (at least in this particular example)
    ! Compare every time?
    ! Some validity check would be nice
    !---------------------------------------------

    class(t_fieldset), intent(in)          :: this
    type(t_datetime), intent(in)           :: date
    character(len=LEN_CHAR_L), intent(out) :: thedir
    integer, optional, intent(out)         :: folder_idx
    integer                                :: i
    integer(rk)                            :: YYYYMMDD
    character(len=LEN_CHAR_L)              :: testdir
#ifdef DEBUG
    integer, save                          :: NCall = 0
#endif

#ifdef DEBUG
    NCall = NCall + 1

    call date%print_short_date()
#endif

    if (this%nentries == 1) then
      write (testdir, '(a,i0.8)') trim(this%PATH)//'/', this%dirlist(1)

      if (date < datetime_from_netcdf(trim(testdir)//PROC0, 1)) then
        call throw_error("fieldset :: find_folder", "The date is earlier than first date in data file")
      end if
      thedir = testdir
      if (present(folder_idx)) folder_idx = 1

      return
    end if

    YYYYMMDD = date%shortDate(.false.)
    do i = 1, this%nentries
      if (this%dirlist(i) .ge. YYYYMMDD) then

        write (testdir, '(a,i0.8)') trim(this%PATH)//'/', this%dirlist(i)

        if (date < datetime_from_netcdf(trim(testdir)//PROC0, 1)) then

          if (i == 1) call throw_error("fieldset :: find_folder", "The date is earlier than first date in data file")
          write (thedir, '(a,i0.8)') trim(this%PATH)//'/', this%dirlist(i - 1)
          if (present(folder_idx)) folder_idx = i - 1

        else

          write (thedir, '(a,i0.8)') trim(this%PATH)//'/', this%dirlist(i)
          if (present(folder_idx)) folder_idx = i

        end if

        return
      end if
    end do

    call throw_error("fieldset :: find_folder", "Did not find folder")

    return
  end subroutine find_folder
  !===========================================
  subroutine find_file(this, date, thefile, file_idx)
    !---------------------------------------------
    ! Assumes a sorted filelist.
    !---------------------------------------------

    class(t_fieldset), intent(in)          :: this
    type(t_datetime), intent(in)           :: date
    character(len=LEN_CHAR_L), intent(out) :: thefile
    integer, optional, intent(out)         :: file_idx
    integer                                :: i, n_times
#ifdef DEBUG
    integer, save                          :: NCall = 0
    type(t_datetime)                       :: testdate
#endif

#ifdef DEBUG
    NCall = NCall + 1

    call date%print_short_date()
#endif

    if (this%nentries == 1) then
      write (thefile, '(a)') trim(this%PATH)//'/'//trim(this%filelist(1))
#ifdef DEBUG
      testdate = datetime_from_netcdf(trim(thefile), 1)

#endif
      if ((date < datetime_from_netcdf(trim(thefile), 1))) then
        call throw_error("fieldset :: find_file", "The date is earlier than first date in data file")
      end if
      if (present(file_idx)) file_idx = 1

      return
    end if

    do i = 1, this%nentries
      if (date < datetime_from_netcdf(trim(this%PATH)//'/'//trim(this%filelist(i)), 1)) then
        if (i == 1) call throw_error("fieldset :: find_file", "The date is earlier than first date in data file")
        write (thefile, '(a)') trim(this%PATH)//'/'//trim(this%filelist(i - 1))
        if (present(file_idx)) file_idx = i - 1

        return
      end if
    end do
    ! If we reach here, just return the last file and hope for the best
    write (thefile, '(a)') trim(this%PATH)//'/'//trim(this%filelist(this%nentries))
    if (date > datetime_from_netcdf(trim(thefile), n_times)) call throw_error("fieldset :: find_file", "Did not find file")
    if (present(file_idx)) file_idx = this%nentries

    return
  end subroutine find_file
  !===========================================
  subroutine get_indices_vertical(this, t, z, i, j, k, kr)
    class(t_fieldset), intent(in)   :: this
    real(rk), intent(in)            :: t, z
    integer, intent(in)             :: i, j
    integer, optional, intent(out)  :: k
    real(rk), optional, intent(out) :: kr
    real(rk)                        :: zax(this%nz)
    real(rk)                        :: dep
    integer                         :: ik

    dep = this%domain%get_bathymetry(i, j)

    if (dep < ZERO) then
      if (present(k)) k = this%nz
      if (present(kr)) kr = real(this%nz, rk)
#ifdef DEBUG
      call throw_warning("fieldset :: get_indices_vertical", "Particle on ground") ! Better error message
      debug(dep)
#endif
      return
    end if

    if (z < -1.0 * dep) then
      if (present(k)) k = 1
      if (present(kr)) kr = 1.0d0
#ifdef DEBUG
      call throw_warning("fieldset :: get_indices_vertical", "Out of bounds") ! Better error message
      debug(dep)
#endif
      return
    end if

    zax = this%get_zax(t, i, j) ! Don't want to interpolate the Z axis in space, so we're taking the closest indices
    
    if (z > zax(this%nz)) then
      if (present(k)) k = this%nz
      if (present(kr)) kr = real(this%nz, kind=rk)
#ifdef DEBUG
      call throw_warning("fieldset :: get_indices_vertical", "Out of bounds") ! Better error message
      debug(dep)
#endif
      return
    end if

    !---------------------------------------------
    ! Could there be a way to do this without looping?
    do ik = 1, this%nz - 1
      if (zax(ik + 1) .ge. z) then
        if ((zax(ik + 1) - zax(ik)) <= ZERO) cycle ! There may be zeros at the bottom
        if (present(k)) then
          k = ik
        end if
        if (present(kr)) then
          kr = ik + (z - zax(ik)) / (zax(ik + 1) - zax(ik))
          if (kr < ONE) kr = ONE
        end if
        return
      end if
    end do
    if (present(k)) then
      k = this%nz
    end if
    if (present(kr)) then
      kr = real(this%nz, rk)
    end if

#ifdef DEBUG
    call throw_warning("fieldset :: get_indices_vertical", "Did not find vertical index!")
#endif

    return
  end subroutine get_indices_vertical
  !===========================================
  subroutine search_indices(this, t, x, y, z, i, j, k, ir, jr, kr)
    class(t_fieldset), intent(in)   :: this
    real(rk), optional, intent(in)  :: t           ! time
    real(rk), intent(in)            :: x, y        ! and position (lon, lat, depth)
    real(rk), optional, intent(in)  :: z
    integer, optional, intent(out)  :: i, j, k     ! integer indices out
    real(rk), optional, intent(out) :: ir, jr, kr  ! real indices out
    integer                         :: ii, jj, kk
    real(rk)                        :: iir, jjr, kkr

    call this%domain%get_indices_2d(x, y, i=ii, j=jj, ir=iir, jr=jjr)
    if (present(k) .or. present(kr)) then
      if (.not. (present(t) .and. present(z))) call throw_error("fieldset :: search_indices", "Time or depth missing!")
      if (this%zax_idx < 1) then ! Probably the fastest way to make sure the Z axis exists
        kk = this%nz; kkr = real(this%nz, kind=rk)
      else
        call this%get_indices_vertical(t, z, ii, jj, k=kk, kr=kkr)
      end if
    end if
    if (present(i)) then
      i = ii; 
    end if
    if (present(j)) then
      j = jj; 
    end if
    if (present(k)) then
      k = kk; 
    end if
    if (present(ir)) then
      ir = iir; 
    end if
    if (present(jr)) then
      jr = jjr; 
    end if
    if (present(kr)) then
      kr = kkr; 
    end if

    return
  end subroutine search_indices
  !===========================================
  function get_zax(this, t, i, j) result(res)
    class(t_fieldset), intent(in) :: this
    real(rk), intent(in)          :: t
    integer, intent(in)           :: i, j
    type(t_field), pointer        :: p_field
    real(rk)                      :: arr_zax(this%nz), tmp_zax(this%nz)
    real(rk)                      :: res(this%nz)
    integer                       :: ik

    call this%fields%get_item(this%zax_idx, p_field)
    call p_field%get_array_1D(t=t, i=i, j=j, n=this%nz, dim=3, res=arr_zax)

    select case (this%zax_style)
    case (DEPTH_VALUES)

      res = arr_zax
    case (LAYER_THICKNESS)

      tmp_zax = arr_zax
      where (tmp_zax <= -5.0d0) tmp_zax = 0.0
      tmp_zax(1) = -1.0 * this%domain%get_bathymetry(i, j)
      do ik = 2, this%nz
        tmp_zax(ik) = tmp_zax(ik - 1) + arr_zax(ik)
      end do
      res = tmp_zax
    case default
      call throw_error("fieldset :: get_zax", "Z - axis style unknown: zax_style = 1 or 2")
    end select

    return
  end function get_zax
  !===========================================
  function get_gradient(this, field_name, t, dim) result(res)
    class(t_fieldset), intent(in) :: this
    character(*), intent(in)      :: field_name
    real(rk), intent(in)          :: t
    integer, intent(in)           :: dim
    type(t_field), pointer        :: p_field
    real(rk)                      :: res(this%nx, this%ny, this%nz)

    call this%fields%get_item(field_name, p_field)
    select case (dim)
    case (1)
      res = p_field%gradient(t, this%domain%dx, dim)
    case (2)
      res = p_field%gradient(t, this%domain%dy, dim)
    case default
      call throw_error("fieldset :: get_gradient", "Dimension must be 1 or 2")
    end select

  end function get_gradient
  !===========================================
  real(rk) function sealevel(this, t, i, j) result(res)
    !---------------------------------------------
    ! Should I keep this here?
    !---------------------------------------------

    class(t_fieldset), intent(in) :: this
    real(rk), intent(in)          :: t, i, j
    type(t_field), pointer        :: p_field
    real(rk)                      :: zax(this%nz)

    res = 0.
    if (this%fields%key_exists("ELEV")) then
      call this%fields%get_item("ELEV", p_field)
      call p_field%get(t, i, j, &
                       res=res)
    else
      zax = this%get_zax(t, nint(i), nint(j))
      res = zax(this%nz)
    end if

    return
  end function sealevel
  !===========================================
  subroutine update(this, date, ignore_check, update_dates)
    class(t_fieldset), intent(inout) :: this
    class(t_datetime), intent(in)    :: date
    logical, optional, intent(in)    :: ignore_check, update_dates
    type(t_field), pointer           :: p_field
    integer                          :: i_field
    logical                          :: ign_chk, ud

    if (present(ignore_check)) then
      ign_chk = ignore_check
    else
      ign_chk = .false.
    end if

    ! Check if it's even time to read
    if ((.not. ign_chk) .and. (date < this%next_read_dt)) then

      return
    end if

#ifdef DEBUG

    call date%print_short_date()

#endif

    ! Read nc (if - else because loop over subdomains might need to be the other way)
    if (this%has_subdomains) then
      do i_field = 1, this%num_fields
        call this%fields%get_item(i_field, p_field)
        call p_field%swap()
        call this%read_field_subdomains(p_field, i_field)
      end do
    else
      do i_field = 1, this%num_fields
        call this%fields%get_item(i_field, p_field)
        call p_field%swap()
        call this%read_field(p_field, i_field)
      end do
    end if

    this%read_idx = this%read_idx + this%read_idx_increment
    if (this%read_idx > this%current_ntimes) then
      this%read_idx = this%read_idx - this%current_ntimes
      ! Hopefully noone will be skipping whole files
      this%dirlist_idx = this%dirlist_idx + 1
      if (this%has_subdomains) then
        this%current_path = this%get_folder(this%dirlist_idx)
        call nc_get_dim(trim(this%current_path)//PROC0, 'time', this%current_ntimes)
      else
        this%current_path = this%get_file(this%dirlist_idx)
        call nc_get_dim(trim(this%current_path), 'time', this%current_ntimes)
      end if
    end if

    if (present(update_dates)) then
      ud = update_dates
    else
      ud = .true.
    end if

    if (.not. ign_chk) call this%next_read_dt%update(this%nc_timestep)
    if (ud) then
      call this%date_t1%update(this%nc_timestep)
      call this%date_t2%update(this%nc_timestep)
    end if

    return
  end subroutine update
  !===========================================
  subroutine read_first_timesteps(this, date)
    class(t_fieldset), intent(inout) :: this
    type(t_datetime), intent(in)     :: date

    call this%update(date, update_dates=.false.)
    call this%update(date, ignore_check=.true., update_dates=.false.)

    return
  end subroutine read_first_timesteps
  !===========================================
  subroutine read_field_subdomains(this, p_field, field_idx)
    class(t_fieldset), intent(in)                :: this
    type(t_field), target, intent(inout)         :: p_field
    integer, intent(in)                          :: field_idx
    character                                    :: c_field
    character(len=LEN_CHAR_S)                    :: varname
    character(len=LEN_CHAR_L)                    :: subdom_filename
    integer                                      :: n_dims
    real(rk), dimension(:, :, :), allocatable    :: buffer, tmp_arr
    integer, allocatable                         :: start(:), count(:)
    integer                                      :: i, j, i_subdom, ioff, joff, istart, jstart

    if (field_idx == this%u_idx) then
      c_field = "u"
    else if (field_idx == this%v_idx) then
      c_field = "v"
    else
      c_field = "x"
    end if

    varname = p_field%get_varname()
    n_dims = p_field%n_dims

    if (n_dims == 2) then
      allocate (buffer(this%nx, this%ny, 1))
      allocate (count(3))
      count(3) = 1
      start = [1, 1, this%read_idx]
    else if (n_dims == 3) then
      allocate (buffer(this%nx, this%ny, this%nz))
      allocate (count(4))
      count(3) = this%nz
      count(4) = 1
      start = [1, 1, 1, this%read_idx]
    else
      call throw_error("fieldset :: read_field", "Wrong number of dimensions: "//trim(varname))
    end if

    ! Makes the uninitialised ubound warning go away, but might not be necessary actually
    buffer = ZERO

    do i_subdom = 0, this%nproc - 1
      write (subdom_filename, "(a,i0.4,a)") trim(this%current_path)//"/"//trim(this%file_prefix), &
        i_subdom, trim(this%file_suffix)//".nc"

      ioff = this%pmap(i_subdom + 1, 1); joff = this%pmap(i_subdom + 1, 2)
      istart = 1 + ioff; jstart = 1 + joff
      ! count = [this%nxp, this%nyp, this%nz, 1]
      if (ioff .lt. 0) then
        count(1) = this%nxp - abs(ioff)
        istart = 1
      else
        count(1) = this%nxp
      end if
      if (joff .lt. 0) then
        count(2) = this%nyp - abs(joff)
        jstart = 1
      else
        count(2) = this%nyp
      end if

#ifdef DEBUG
      if (mod(i_subdom, 50) .eq. 0) then

      end if
#endif

      allocate (tmp_arr(count(1), count(2), count(3)))

      if (n_dims == 2) then
        call nc_read_real_3d(trim(subdom_filename), trim(varname), start, count, tmp_arr)
      else
        call nc_read_real_4d(trim(subdom_filename), trim(varname), start, count, tmp_arr)
      end if
      buffer(istart:istart + count(1) - 1, jstart:jstart + count(2) - 1, 1:count(3)) = tmp_arr

      deallocate (tmp_arr)

    end do

    select case (c_field)
    case ("v")
      do j = 1, this%ny
        do i = 1, this%nx
          if (this%domain%get_bathymetry(i, j) < ZERO) then
            buffer(i, j, :) = ZERO
          end if
          ! This is needed due to the A grid
          if (this%u_mask(i, j) > 0) then
            buffer(i, j, :) = ZERO
          end if
        end do
      end do
      do j = 1, this%ny
        do i = 1, this%nx
          if (this%v_mask(i, j) > 0) then
            buffer(i, j, :) = buffer(i + 1, j, :)
          else if (this%v_mask(i, j) < 0) then
            buffer(i + 1, j, :) = buffer(i, j, :)
          end if
        end do
      end do
    case ("u")
      do j = 1, this%ny
        do i = 1, this%nx
          if (this%domain%get_bathymetry(i, j) < ZERO) then
            buffer(i, j, :) = ZERO
          end if
          if (this%v_mask(i, j) > 0) then
            buffer(i, j, :) = ZERO
          end if
        end do
      end do
      do j = 1, this%ny
        do i = 1, this%nx
          if (this%u_mask(i, j) > 0) then
            buffer(i, j, :) = buffer(i, j + 1, :)
          else if (this%u_mask(i, j) < 0) then
            buffer(i, j + 1, :) = buffer(i, j, :)
          end if
        end do
      end do
    case default
      do j = 1, this%ny
        do i = 1, this%nx
          if (this%domain%get_bathymetry(i, j) < ZERO) then
            buffer(i, j, :) = ZERO
          end if
        end do
      end do
    end select

    where (buffer <= MISSING_VAL) buffer = ZERO

    call p_field%set(buffer)

  end subroutine read_field_subdomains
  !===========================================
  subroutine read_field(this, p_field, field_idx)
    class(t_fieldset), intent(in)             :: this
    type(t_field), pointer, intent(inout)     :: p_field
    integer, intent(in)                       :: field_idx
    character                                 :: c_field
    character(len=LEN_CHAR_S)                 :: varname
    integer                                   :: n_dims
    real(rk), dimension(:, :, :), allocatable :: buffer
    integer, allocatable                      :: start(:), count(:)
    integer                                   :: i, j

    if (field_idx == this%u_idx) then
      c_field = "u"
    else if (field_idx == this%v_idx) then
      c_field = "v"
    else
      c_field = "x"
    end if

    varname = p_field%get_varname()
    n_dims = p_field%n_dims

    if (n_dims == 2) then
      start = [1, 1, this%read_idx]
      count = [this%nx, this%ny, 1]
    else if (n_dims == 3) then
      start = [1, 1, 1, this%read_idx]
      count = [this%nx, this%ny, this%nz, 1]
    else
      call throw_error("fieldset :: read_field", "Wrong number of dimensions: "//trim(varname))
    end if

    allocate (buffer(count(1), count(2), count(3)))

    if (n_dims == 2) then
      call nc_read_real_3d(trim(this%current_path), trim(varname), start, count, buffer)
    else
      call nc_read_real_4d(trim(this%current_path), trim(varname), start, count, buffer)
    end if

    select case (c_field)
    case ("v")
      do j = 1, this%ny
        do i = 1, this%nx
          if (this%domain%get_bathymetry(i, j) < ZERO) then
            buffer(i, j, :) = ZERO
          end if
          ! This is needed due to the A grid
          if (this%u_mask(i, j) > 0) then
            buffer(i, j, :) = ZERO
          end if
        end do
      end do
      do j = 1, this%ny
        do i = 1, this%nx
          if (this%v_mask(i, j) > 0) then
            buffer(i, j, :) = buffer(i + 1, j, :)
          else if (this%v_mask(i, j) < 0) then
            buffer(i + 1, j, :) = buffer(i, j, :)
          end if
        end do
      end do
    case ("u")
      do j = 1, this%ny
        do i = 1, this%nx
          if (this%domain%get_bathymetry(i, j) < ZERO) then
            buffer(i, j, :) = ZERO
          end if
          if (this%v_mask(i, j) > 0) then
            buffer(i, j, :) = ZERO
          end if
        end do
      end do
      do j = 1, this%ny
        do i = 1, this%nx
          if (this%u_mask(i, j) > 0) then
            buffer(i, j, :) = buffer(i, j + 1, :)
          else if (this%u_mask(i, j) < 0) then
            buffer(i, j + 1, :) = buffer(i, j, :)
          end if
        end do
      end do
    case default
      do j = 1, this%ny
        do i = 1, this%nx
          if (this%domain%get_bathymetry(i, j) < ZERO) then
            buffer(i, j, :) = ZERO
          end if
        end do
      end do
    end select

    where (buffer <= MISSING_VAL) buffer = ZERO

    call p_field%set(buffer)

  end subroutine read_field
end module mod_fieldset
