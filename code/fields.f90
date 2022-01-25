#include "cppdefs.h"
module fields
  !----------------------------------------------------------------
  ! This module contains methods for the domain and subdomains.
  ! Domain is initialised in init_domain.
  ! Includes:
  ! - find the folder with the right date
  ! - find the right subdomain (TODO: if subdomains are used)
  !----------------------------------------------------------------
  use precdefs
  use errors
  use params, only: mu, do_velocity
  use domain_vars, only: nx, ny, depdata, x0, y0, dx, dy
  use field_vars
  use modtime
  use time_vars, only: run_start_dt
  use nc_manager, only: nc_read4d, nc_var_exists
  use physics, only: seawater_density_from_temp_and_salt
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: init_dirlist, return_pmap, init_fields, find_folder, find_file, &
            get_indices2d, get_indices_vert, find_proc, &
            read_fields_full_domain, &
            sealevel
  !---------------------------------------------
  ! Directories
  character(len=*), parameter     :: dirinfile = 'dirlist.dat' ! List the dirs into this file
  character(len=128), allocatable :: filelist(:)               ! List of files
  integer, allocatable            :: dirlist(:)                ! List of directories (one or the other will be allocated)
  integer                         :: nentries                  ! Number of directories or files
  !---------------------------------------------
  ! Subdomain variables
  integer              :: nproc      ! Number of subdomains
  integer              :: nxp, nyp   ! Size of subdomain
  integer, allocatable :: pmap(:, :)  ! Offsets of subdomains
  integer, allocatable :: pmask(:, :) ! Proc. mask
  !---------------------------------------------
  integer :: ierr
  !===================================================
contains
  !===========================================
  subroutine init_dirlist
    !---------------------------------------------
    ! Get a list of directories or files that contain the data
    ! The files should be named so that ls command would give them
    ! in the right order. Including the date in the file name should be enough.
    ! TODO: Alternatively, if all the files are in one directory,
    !       it should list all the files in GETMPATH.
    !       (separate routine e.g., init_filelist ?)
    ! TODO (later, probably never): sort the files somehow so all of this
    !       would not depend on ls getting it right.
    !---------------------------------------------

    logical :: dirlist_exists
    integer :: idir

    FMT1, "======== Init dirlist ========"
    inquire (file=trim(dirinfile), exist=dirlist_exists)
    if (.not. dirlist_exists) then
#ifndef NOSYSCALLS
      select case (has_subdomains)
      case (.true.)
        FMT2, "Getting directory list..."
        call system('( ls -d '//trim(GETMPATH)//'/*/ | wc -l && ls -d ' &
                    //trim(GETMPATH)//'/*/ | xargs -n 1 basename ) > ' &
                    //trim(dirinfile), status=ierr)
        if (ierr .ne. 0) call throw_error("init_dirlist", "Could not get directory list from "//trim(GETMPATH), ierr)
      case (.false.)
        FMT2, "Getting file list..."
        call system('( ls '//trim(GETMPATH)//'/'//trim(file_prefix)//'*' &
                    //trim(file_suffix)//'*.nc | wc -l && ls '//trim(GETMPATH) &
                    //'/'//trim(file_prefix)//'*'//trim(file_suffix) &
                    //'*.nc | xargs -n 1 basename ) > '//trim(dirinfile), status=ierr)
        if (ierr .ne. 0) call throw_error("init_dirlist", "Could not get file list from "//trim(GETMPATH), ierr)
      end select
#else
      call throw_error("init_dirlist", "File/directory list ("//trim(dirinfile)//") does not exist!")
#endif
    end if

    open (DIRFILE, file=trim(dirinfile), action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error("init_dirlist", "Failed to open "//trim(dirinfile), ierr)
    read (DIRFILE, *) nentries
    FMT2, "Found ", nentries, "files/directories in "//trim(GETMPATH)
    select case (has_subdomains)
    case (.true.)
      allocate (dirlist(nentries))
      do idir = 1, nentries
        read (DIRFILE, *) dirlist(idir)
      end do
    case (.false.)
      allocate (filelist(nentries))
      do idir = 1, nentries
        read (DIRFILE, '(a)') filelist(idir)
      end do
    end select
    close (DIRFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error("init_dirlist", "Failed to close "//trim(dirinfile), ierr)

    FMT2, "Dirlist initialized"

  end subroutine init_dirlist
  !===========================================
  subroutine get_proc_mask
    !---------------------------------------------
    ! This maps the pieces of GETM data using par_setup.
    ! TODO: Since pmask takes up a lot of unnecessary space,
    !       maybe this should only be called when par_setup is wanted
    !      (either a compilation flag or "if (parallel)" or something in init_model)
    ! EDIT: This is called from init_model only if has_subdomains=.true.
    !---------------------------------------------

    integer :: imax, jmax ! Total size of the domain
    integer :: ioff, joff ! Subdomain offset
    integer :: pnum, iend, jend
    integer :: i, j, iproc

    FMT1, "======== Init subdomains ========"

    open (PROCFILE, file=trim(PMAPFILE), action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error("get_proc_mask", "Failed to open "//trim(PMAPFILE), ierr)
    read (PROCFILE, *) nproc
    allocate (pmap(nproc, 2))
    read (PROCFILE, *) nxp, nyp, imax, jmax
    if ((imax .ne. nx) .or. (jmax .ne. ny)) then
      ERROR, "(nx = ", nx, "imax = ", imax, "ny = ", ny, "jmax = ", jmax, ")"
      call throw_error("get_proc_mask", "Given domain size does not match one in parallel setup file!")
    else
      FMT2, "Number of subdomains: ", nproc, " (", nxp, " x ", nyp, ")"
    end if
    allocate (pmask(imax, jmax))
    do i = 1, imax
      do j = 1, jmax
        pmask(i, j) = -10
      end do
    end do
    iproc = 0
    do while (iproc .lt. nproc)
      read (PROCFILE, *) pnum, ioff, joff
      i = ioff + 1; j = joff + 1
      iproc = iproc + 1
      pmap(iproc, 1) = ioff; pmap(iproc, 2) = joff
      iend = nxp; jend = nyp
      if (i .lt. 1) then
        !---------------------------------------------
        ! If i is less than 1, then ioff must be negative.
        ! Therefore this processor actually covers less area,
        ! i.e. iend must be reduced ioff amount
        iend = nxp + ioff ! Same as nxp - abs(ioff)
        i = 1
      end if
      if (j .lt. 1) then
        !---------------------------------------------
        ! Same logic here
        jend = nyp + joff
        j = 1
      end if
      pmask(i:i + iend, j:j + jend) = pnum
    end do
    close (PROCFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error("get_proc_mask", "Failed to close "//trim(PMAPFILE), ierr)

    FMT2, "Processor mask initialised"

    return
  end subroutine get_proc_mask
  !===========================================
  subroutine return_pmap(pmapout)
    !---------------------------------------------
    ! Just in case I want to use pmap somewhere else
    ! (it's private otherwise)
    !---------------------------------------------

    integer, intent(inout) :: pmapout(nproc, 2)

    pmapout = pmap

    return
  end subroutine return_pmap
  !===========================================
  subroutine init_fields
    !---------------------------------------------
    ! Allocate arrays for current data.
    ! TODO: Right now it is assumed that all subdomains
    !       are the same size. Also it is assumed that subdomains
    !       exist at all. There should be a switch (e.g. has_subdomains).
    !       This also changes allocation.
    ! EDIT: Only full domain is used, even when the data is in chunks,
    !       since reading into subdomains is not ready yet.
    !---------------------------------------------
    character(len=256) :: initPath
    character(len=516) :: filename

    FMT1, "======== Init fields ========"

    !---------------------------------------------
    ! Arrays for h. velocity
    FMT2, "Using full domain"
    FMT2, "Allocating fields of size (nx, ny, nz): (", nx, ", ", ny, ", ", nlevels, ")"
    allocate (udata(nx, ny, nlevels), vdata(nx, ny, nlevels), &
              udatanew(nx, ny, nlevels), vdatanew(nx, ny, nlevels), stat=ierr)
    if (ierr .ne. 0) call throw_error("init_fields", "Could not allocate", ierr)
    udata = 0; vdata = 0; udatanew = 0; vdatanew = 0; 
    allocate (udata_interp(nx, ny, nlevels), vdata_interp(nx, ny, nlevels), &
    udatanew_interp(nx, ny, nlevels), vdatanew_interp(nx, ny, nlevels), stat=ierr)
    if (ierr .ne. 0) call throw_error("init_fields", "Could not allocate", ierr)
    udata_interp = 0; vdata_interp = 0; udatanew_interp = 0; vdatanew_interp = 0; 
    !---------------------------------------------
    ! Arrays for v. velocity
    if (run_3d) then
      allocate (wdata(nx, ny, nlevels), wdatanew(nx, ny, nlevels), &
      zaxdata(nx, ny, nlevels), zaxdatanew(nx, ny, nlevels), stat=ierr)
      if (ierr .ne. 0) call throw_error("init_fields", "Could not allocate", ierr)
      wdata = 0; wdatanew = 0; zaxdata = 0; zaxdatanew = 0; 
      allocate (wdata_interp(nx, ny, nlevels), wdatanew_interp(nx, ny, nlevels), &
      zaxdata_interp(nx, ny, nlevels), zaxdatanew_interp(nx, ny, nlevels), stat=ierr)
      if (ierr .ne. 0) call throw_error("init_fields", "Could not allocate", ierr)
      wdata_interp = 0; wdatanew_interp = 0; zaxdata_interp = 0; zaxdatanew_interp = 0; 
    end if
    !---------------------------------------------
    if (do_velocity) then
      !---------------------------------------------
      ! Arrays for density
      select case (has_subdomains)
      case (.true.)
        call find_folder(run_start_dt, initPath)
        write (filename, '(a)') trim(initPath)//PROC0
      case (.false.)
        call find_file(run_start_dt, initPath)
        write (filename, '(a)') trim(initPath)
      end select
      if (nc_var_exists(trim(filename), "rho")) then
        allocate (density(nx, ny, nlevels), densitynew(nx, ny, nlevels), stat=ierr)
        if (ierr .ne. 0) call throw_error("init_fields", "Could not allocate", ierr)
        density = 0; densitynew = 0
        has_density = DENSITY
      else
        call throw_warning("init_fields", "Could not find density ('rho') in "//trim(filename))
        if (nc_var_exists(trim(filename), "temp") .and. &
        nc_var_exists(trim(filename), "salt")) then
          allocate (temp(nx, ny, nlevels), tempnew(nx, ny, nlevels), stat=ierr)
          if (ierr .ne. 0) call throw_error("init_fields", "Could not allocate", ierr)
          temp = 0; tempnew = 0
          allocate (salt(nx, ny, nlevels), saltnew(nx, ny, nlevels), stat=ierr)
          if (ierr .ne. 0) call throw_error("init_fields", "Could not allocate", ierr)
          salt = 0; saltnew = 0
          has_density = TEMP_SALT
        else
          call throw_warning("init_fields", "Could not find temperature or salinity ('temp'/'salt') in "//trim(filename)// &
          ". Using default density.")
          has_density = DEFAULT_DENSITY
        end if
      end if
      !---------------------------------------------
      ! Arrays for viscosity
      if (nc_var_exists(trim(filename), "nuh")) then
        allocate (visc(nx, ny, nlevels), viscnew(nx, ny, nlevels), stat=ierr)
        if (ierr .ne. 0) call throw_error("init_fields", "Could not allocate", ierr)
        visc = 0; viscnew = 0
        has_viscosity = .true.
      else
        call throw_warning("init_fields", "Could not find viscosity ('nuh') in "//trim(filename)// &
                           ". Using default viscosity.")
      end if
    end if
    !---------------------------------------------
    if (has_subdomains) call get_proc_mask

    FMT2, "Fields allocated"
    !---------------------------------------------
    ! For later...
    ! if (do_subdom) then
    !     FMT2, "Using subdomains"
    !     FMT2, "Allocating fields of size (nx, ny, nz): (", nxp, ", ", nyp, ", ", nlevels, ")"
    !     allocate(udata(nxp, nyp, nlevels), vdata(nxp, nyp, nlevels), &
    !              udatanew(nxp, nyp, nlevels), vdatanew(nxp, nyp, nlevels))
    !     if (run_3d) then
    !       allocate(wdata(nxp, nyp, nlevels), wdatanew(nxp, nyp, nlevels), &
    !                zaxdata(nxp, nyp, nlevels), zaxdatanew(nxp, nyp, nlevels))
    !     end if
    ! end if
    !---------------------------------------------

  end subroutine init_fields
  !===========================================
  subroutine find_folder(date, thedir)
    !---------------------------------------------
    ! TODO: might have to check if nc files start at time
    ! [date] 00:00:00 or [date] 00:00:10 (at least in this particular example)
    ! Compare every time?
    ! Some validity check would be nice
    !---------------------------------------------

    integer                         :: i
    integer, save                   :: NCall = 0
    type(datetime), intent(in)      :: date
    integer(rk)                         :: YYYYMMDD
    character(len=256)              :: testdir
    character(len=256), intent(out) :: thedir

    dbghead(find_folder)

    NCall = NCall + 1
#ifdef DEBUG
    debug(NCall)
    DBG, "Date in:"
    call date%print_short_date()
#endif

    YYYYMMDD = date%shortDate(.false.)

    if (nentries == 1) then
      write (testdir, '(a,i0.8)') trim(GETMPATH)//'/', dirlist(1)
      debug(trim(testdir))
      if (date < init_datetime_from_netcdf(trim(testdir)//PROC0, 1)) then
        call throw_error("find_folder", "The date is earlier than first date in data file")
      end if
      thedir = testdir
      dbgtail(find_folder)
      return
    end if

    do i = 1, nentries
      if (dirlist(i) .ge. YYYYMMDD) then
        debug(YYYYMMDD)
        debug(dirlist(i))
        if (NCall .le. 3) then ! It works, shut up
          DBG, "NCall <= 3"
          write (testdir, '(a,i0.8)') trim(GETMPATH)//'/', dirlist(i)
          debug(testdir)
          if (date < init_datetime_from_netcdf(trim(testdir)//PROC0, 1)) then
            DBG, "date < init_datetime_from_netcdf(trim(testdir)//PROC0, 1)"
            write (thedir, '(a,i0.8)') trim(GETMPATH)//'/', dirlist(i - 1)
            debug(thedir)
          else
            DBG, "date > init_datetime_from_netcdf(trim(testdir)//PROC0, 1)"
            write (thedir, '(a,i0.8)') trim(GETMPATH)//'/', dirlist(i)
            debug(thedir)
          end if
        else ! Probably could be condensed (if and else are exactly the same, but I'm not ready to delete anything)
          DBG, "NCall > 3"
          write (testdir, '(a,i0.8)') trim(GETMPATH)//'/', dirlist(i)
          debug(testdir)
          if (date < init_datetime_from_netcdf(trim(testdir)//PROC0, 1)) then
            DBG, "date < init_datetime_from_netcdf(trim(testdir)//PROC0, 1)"
            write (thedir, '(a,i0.8)') trim(GETMPATH)//'/', dirlist(i - 1)
          else
            DBG, "date > init_datetime_from_netcdf(trim(testdir)//PROC0, 1)"
            write (thedir, '(a,i0.8)') trim(GETMPATH)//'/', dirlist(i)
          end if
          debug(thedir)
        end if
        dbgtail(find_folder)
        return
      end if
    end do

    dbgtail(find_folder)
    return
  end subroutine find_folder
  !===========================================
  subroutine find_file(date, thefile)
    !---------------------------------------------
    ! Assumes a sorted filelist.
    !---------------------------------------------

    integer                         :: i
    integer, save                   :: NCall = 0
    type(datetime), intent(in)      :: date
    type(datetime)                  :: testdate
    character(len=256), intent(out) :: thefile

    dbghead(find_file)

    NCall = NCall + 1
#ifdef DEBUG
    debug(NCall)
    DBG, "Date in:"
    call date%print_short_date()
#endif

    if (nentries == 1) then
      write (thefile, '(a)') trim(GETMPATH)//'/'//trim(filelist(1))
      testdate = init_datetime_from_netcdf(trim(thefile), 1)
#ifdef DEBUG
      DBG, "Test date:"; call testdate%print_short_date()
#endif
      if ((NCall == 1) .and. (date < testdate)) then
        call throw_error("find_file", "The date is earlier than first date in data file")
      end if
      dbgtail(find_file)
      return
    end if

    do i = 1, nentries
      if (date < init_datetime_from_netcdf(trim(GETMPATH)//'/'//trim(filelist(i)), 1)) then
        write (thefile, '(a)') trim(GETMPATH)//'/'//trim(filelist(i - 1))
        dbgtail(find_file)
        return
      end if
    end do
    ! If we reach here, just return the last file and hope for the best
    write (thefile, '(a)') trim(GETMPATH)//'/'//trim(filelist(nentries))

    dbgtail(find_file)
    return
  end subroutine find_file
  !===========================================
  subroutine get_indices2d(xin, yin, x0in, y0in, dxin, dyin, i, j, ir, jr)

    integer, intent(out)            :: i, j
    real(rk), intent(in)            :: xin, yin, x0in, y0in, dxin, dyin
    real(rk), intent(out), optional :: ir, jr
    real(rk)                        :: irt, jrt

    irt = (xin - x0in) / dxin; jrt = (yin - y0in) / dyin

    i = int(irt); j = int(jrt)

    ! TODO: Throw an error instead?
    if (i .lt. 1) i = 1
    if (j .lt. 1) j = 1

    if (present(ir) .and. present(jr)) then
      ir = irt; jr = jrt
    end if

    return
  end subroutine get_indices2d
  !===========================================
  !subroutine get_indices_vert(zin, z0in, dzin, k, kr)
  subroutine get_indices_vert(zaxarr, zin, i, j, k, kr, dzout)
    !---------------------------------------------
    ! Get 3D indices.
    ! Uses "zcn" from GETM output and assumes constant dz.
    ! TODO: Maybe it would be better to use layer thickness instead? (later)
    ! Went with another approach, this now also returns dz...
    !---------------------------------------------

    integer, intent(in)             :: i, j
    integer, intent(out)            :: k
    integer                         :: ik
    real(rk), intent(in)            :: zin
    real(rk), intent(in)            :: zaxarr(nx, ny, nlevels)
    real(rk), intent(out), optional :: kr, dzout
    real(rk)                        :: tmp_zax(nlevels)

    dbghead(get_indices_vert)

    debug(zin)

    select case (zax_style)
    case (DEPTH_VALUES)
      DBG, "Using depth values"
      tmp_zax = zaxarr(i, j, :)
    case (LAYER_THICKNESS)
      DBG, "Using layer thickness"
      tmp_zax = zaxarr(i, j, :)
      debug(tmp_zax)
      where (tmp_zax <= -5.0d0) tmp_zax = 0.0
      tmp_zax(1) = -1.0 * depdata(i, j)
      do ik = 2, nlevels
        tmp_zax(ik) = tmp_zax(ik - 1) + zaxarr(i, j, ik)
      end do
    case default
      call throw_error("get_indices_vert", "Z - axis style unknown: zax_style = 1 or 2")
    end select
    debug(tmp_zax)

    do ik = 1, nlevels - 1
      if (tmp_zax(ik + 1) .ge. zin) then
        k = ik
        if (present(kr)) then
          kr = k + (zin - tmp_zax(ik)) / (tmp_zax(ik + 1) - tmp_zax(ik))
        end if
        if (present(dzout)) then
          dzout = abs(tmp_zax(ik) - tmp_zax(ik + 1))
        end if
        if (kr .lt. 1.0d0) kr = 1.0d0
        debug(tmp_zax(ik)); 
        debug(tmp_zax(ik + 1))
        debug(k); debug(kr); debug(dzout)
        dbgtail(get_indices_vert)
        return
      end if
    end do
    call throw_warning("get_indices_vert", "Did not find vertical index!")

    dbgtail(get_indices_vert)
    return
  end subroutine get_indices_vert
  !===========================================
  real(rk) function sealevel(zaxarr, i, j) result(res)

    integer, intent(in)  :: i, j
    integer              :: ik
    real(rk), intent(in) :: zaxarr(nx, ny, nlevels)
    real(rk)             :: tmp_zax(nlevels)

    dbghead(sealevel)
    debug(shape(zaxarr))

    select case (zax_style)
    case (DEPTH_VALUES)
      DBG, "Using depth values"
      tmp_zax = zaxarr(i, j, :)
    case (LAYER_THICKNESS)
      DBG, "Using layer thickness"
      tmp_zax = zaxarr(i, j, :)
      where (tmp_zax <= -5.0d0) tmp_zax = 0.0
      tmp_zax(1) = -1.0 * depdata(i, j)
      do ik = 2, nlevels
        tmp_zax(ik) = tmp_zax(ik - 1) + zaxarr(i, j, ik)
      end do
      debug(tmp_zax)
    case default
      call throw_error("sealevel", "Z - axis style unknown: zax_style = 1 or 2")
    end select
    res = tmp_zax(nlevels)
    debug(res)

    dbgtail(sealevel)
    return
  end function sealevel
  !===========================================
  subroutine find_proc(iglobal, jglobal, ilocal, jlocal, pnum)
    !---------------------------------------------
    ! Get the processor number and local indices
    ! TODO: Error check for pnum = -10
    !---------------------------------------------

    integer, intent(in)  :: iglobal, jglobal
    integer, intent(out) :: ilocal, jlocal, pnum

    pnum = pmask(iglobal, jglobal)
    if (pnum .eq. -10) then
      call throw_warning("find_proc", "Particle not in any subdomain!")
    end if
    ilocal = iglobal - pmap(pnum + 1, 1)
    jlocal = jglobal - pmap(pnum + 1, 2)

    return
  end subroutine find_proc
  !===========================================
  subroutine read_fields_full_domain(timeindex, path, read_first)
    !---------------------------------------------
    ! Just thinking out loud here... Since I want
    ! the velocity field to be the same size as
    ! topo and seamask, and all the subdomains are
    ! the same size (nxp x nyp) and ioff/joff may
    ! be negative, then the offsets in pmap must
    ! be taken into account, i.e. I will read the
    ! whole subdomain and write it into
    ! udata(ioff:ioff+nxp, ...), EXCEPT when the
    ! offset is negative. In that case I can't read
    ! the whole subdomain. I must read
    ! velx(1+abs(ioff):nxp, ...) and write it as
    ! udata(1:nxp-abs(ioff), ...). It should then
    ! cover less area, like pmask.
    ! EDIT: TURNS OUT THE SUBDOMS ARE NOT THE SAME SIZE, WHATTHEHELL
    ! TODO: How to eliminate the warning?
    ! EDIT: The warning disappeared when using buffers,
    !       but is seems to be much slower this way...
    !---------------------------------------------

    logical                               :: vertical_read = .false., &
                                             density_read = .false., &
                                             temp_salt_read = .false., &
                                             viscosity_read = .false.
    logical, intent(in)                   :: read_first
    integer                               :: iproc, &
                                             ioff, joff, &
                                             istart, jstart, &
                                             start(4), count(4), &
                                             ii, jj
    integer, intent(in)                   :: timeindex
    character(len=256), intent(in)        :: path
    character(len=516)                    :: subdom_filename
    real(rk), dimension(:, :, :), pointer :: arr_u, arr_v, &
                                             arr_w, &
                                             arr_zax, &
                                             arr_density, &
                                             arr_temp, &
                                             arr_salt, &
                                             arr_viscosity
    real(rk), allocatable                 :: buffer(:, :, :, :)

    dbghead(read_fields_full_domain)
    debug(has_subdomains)

    vertical_read = run_3d
    viscosity_read = has_viscosity
    select case (has_density)
    case (DENSITY)
      density_read = .true.
    case (TEMP_SALT)
      temp_salt_read = .true.
    end select

    if (.not. read_first) then
      arr_u => udata
      arr_v => vdata
      arr_w => wdata
      arr_zax => zaxdata
      arr_density => density
      arr_temp => temp
      arr_salt => salt
      arr_viscosity => visc
    else
      arr_u => udatanew
      arr_v => vdatanew
      arr_w => wdatanew
      arr_zax => zaxdatanew
      arr_density => densitynew
      arr_temp => tempnew
      arr_salt => saltnew
      arr_viscosity => viscnew
    end if

    select case (has_subdomains)
    case (.true.)
      do iproc = 0, nproc - 1
        write (subdom_filename, "(a,i0.4,a)") trim(path)//"/"//trim(file_prefix), iproc, trim(file_suffix)//".nc"
        ioff = pmap(iproc + 1, 1); joff = pmap(iproc + 1, 2)
        istart = 1 + ioff; jstart = 1 + joff
        start = [1, 1, startlevel, timeindex]
        count = [nxp, nyp, nlevels, 1]
        if (ioff .lt. 0) then
          count(1) = nxp - abs(ioff)
          istart = 1
        end if
        if (joff .lt. 0) then
          count(2) = nyp - abs(joff)
          jstart = 1
        end if

#ifdef DEBUG
        if (mod(iproc, 50) .eq. 0) then
          debug(trim(subdom_filename))
          debug(start)
          debug(count)
        end if
#endif

        allocate (buffer(count(1), count(2), count(3), count(4)))

        call nc_read4d(trim(subdom_filename), trim(uvarname), start, count, buffer)
        arr_u(istart:istart + count(1) - 1, jstart:jstart + count(2) - 1, 1:nlevels) = buffer(:, :, :, 1)
        call nc_read4d(trim(subdom_filename), trim(vvarname), start, count, buffer)
        arr_v(istart:istart + count(1) - 1, jstart:jstart + count(2) - 1, 1:nlevels) = buffer(:, :, :, 1)

        if (vertical_read) then
          call nc_read4d(trim(subdom_filename), trim(wvarname), start, count, buffer)
          arr_w(istart:istart + count(1) - 1, jstart:jstart + count(2) - 1, 1:nlevels) = buffer(:, :, :, 1)
          call nc_read4d(trim(subdom_filename), trim(zaxvarname), start, count, buffer)
          arr_zax(istart:istart + count(1) - 1, jstart:jstart + count(2) - 1, 1:nlevels) = buffer(:, :, :, 1)
        end if

        if (density_read) then
          call nc_read4d(trim(subdom_filename), "rho", start, count, buffer)
          arr_density(istart:istart + count(1) - 1, jstart:jstart + count(2) - 1, 1:nlevels) = buffer(:, :, :, 1)
        end if

        if (temp_salt_read) then
          call nc_read4d(trim(subdom_filename), "temp", start, count, buffer)
          arr_temp(istart:istart + count(1) - 1, jstart:jstart + count(2) - 1, 1:nlevels) = buffer(:, :, :, 1)
          call nc_read4d(trim(subdom_filename), "salt", start, count, buffer)
          arr_salt(istart:istart + count(1) - 1, jstart:jstart + count(2) - 1, 1:nlevels) = buffer(:, :, :, 1)
        end if

        if (viscosity_read) then
          call nc_read4d(trim(subdom_filename), "nuh", start, count, buffer)
          arr_viscosity(istart:istart + count(1) - 1, jstart:jstart + count(2) - 1, 1:nlevels) = buffer(:, :, :, 1)
        end if

        deallocate (buffer)
      end do
    case (.false.)
      write (subdom_filename, '(a)') trim(path)
      start = [1, 1, startlevel, timeindex]
      count = [nx, ny, nlevels, 1]

      debug(start)
      debug(count)

      call nc_read4d(trim(subdom_filename), trim(uvarname), start, count, arr_u)
      call nc_read4d(trim(subdom_filename), trim(vvarname), start, count, arr_v)

      if (vertical_read) then
        call nc_read4d(trim(subdom_filename), trim(wvarname), start, count, arr_w)
        call nc_read4d(trim(subdom_filename), trim(zaxvarname), start, count, arr_zax)
      end if

      if (density_read) then
        call nc_read4d(trim(subdom_filename), "rho", start, count, arr_density)
      end if

      if (temp_salt_read) then
        call nc_read4d(trim(subdom_filename), "temp", start, count, arr_temp)
        call nc_read4d(trim(subdom_filename), "salt", start, count, arr_salt)
      end if

      if (viscosity_read) then
        call nc_read4d(trim(subdom_filename), "nuh", start, count, arr_viscosity)
      end if

    end select

    do ii = 1, nx
      do jj = 1, ny
        if (depdata(ii, jj) .lt. 0.0d0) then
          arr_u(ii, jj, :) = 0.0d0
          arr_v(ii, jj, :) = 0.0d0

          if (vertical_read) then
            arr_w(ii, jj, :) = 0.0d0
            arr_zax(ii, jj, :) = 0.0d0
          end if

          if (density_read) then
            arr_density(ii, jj, :) = 0.0d0
          end if

          if (temp_salt_read) then
            arr_temp(ii, jj, :) = 0.0d0
            arr_salt(ii, jj, :) = 0.0d0
          end if

          if (viscosity_read) then
            arr_viscosity(ii, jj, :) = 0.0d0
          end if

        end if
      end do
    end do

    where (arr_u <= MISSING_VAL) arr_u = 0.0d0
    where (arr_v <= MISSING_VAL) arr_v = 0.0d0
    if (vertical_read) then
      where (arr_w <= MISSING_VAL) arr_w = 0.0d0
      where (arr_zax <= MISSING_VAL) arr_zax = 0.0d0
    end if
    if (density_read) where (arr_density <= MISSING_VAL) arr_density = 0.0d0
    if (temp_salt_read) then
      where (arr_temp <= MISSING_VAL) arr_temp = 0.0d0
      where (arr_salt <= MISSING_VAL) arr_salt = 0.0d0
    end if
    if (viscosity_read) where (arr_viscosity <= MISSING_VAL) arr_viscosity = 0.0d0

    dbgtail(read_fields_full_domain)
    return
  end subroutine read_fields_full_domain

end module fields
