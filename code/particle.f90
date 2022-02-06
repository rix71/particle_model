#include "cppdefs.h"
module particle_type
  !----------------------------------------------------------------
  ! This is the particle type definition
  !----------------------------------------------------------------
  use precdefs
  use domain_vars, only: x0, y0, dx, dy, seamask, depdata, nx, ny
  ! Pass loop vars into functions rather than import?
  use loop_vars, only: ig, jg, igr, jgr, xnew, ynew, znew, &
                       pvelu, pvelunew, pvelv, pvelvnew, pvelw, pvelwnew
  use time_vars, only: dt
  use field_vars, only: nlevels, uspeed, vspeed, wspeed, &
                        uspeednew, vspeednew, wspeednew, &
                        run_3d
  use fields, only: get_indices2d, sealevel
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_particle
  !---------------------------------------------
  ! Particle type
  type t_particle
    logical           :: isActive = .true.     ! Skip particle in loop if isActive == .false.
    logical           :: kill_bch, kill_bdy    ! Set isActive=.false. if beached or on boundary?
    integer           :: warnings = 0
    real              :: state = 0.            ! 0 - active, 1 - beached, 2 - on boundary, 3 - bottom
    character(len=64) :: originName = 'NoName' ! Origin of particle, e.g. city/country name
    real(rk)          :: xPos = 0.0d0          ! Particle position
    real(rk)          :: yPos = 0.0d0          ! Particle position
    real(rk)          :: zPos = 0.0d0          ! Particle position
    real(rk)          :: u = 0.0d0             ! Particle velocity
    real(rk)          :: v = 0.0d0             ! Particle velocity
    real(rk)          :: w = 0.0d0             ! Particle velocity
    real(rk)          :: rho = 0.0d0           ! Particle density
    real(rk)          :: radius = 0.0d0        ! Particle radius
    real(rk)          :: age = 0.0d0           ! Particle age
    real(rk)          :: trajLen = 0.0d0       ! Particle trajectory length
    real(rk)          :: timeOnBeach = 0.0d0   ! Time spent in the beach area
    real(rk)          :: beachingtime          ! Different particles may essentialy have different beaching times
    real(rk)          :: originNum             ! Origin of particle, number

  contains
    procedure :: print_particle_info
    procedure :: update
    procedure :: check_beached_bdy
    procedure :: check_age
    procedure :: check_depth
  end type t_particle

  !===================================================
contains
  !===========================================
  subroutine print_particle_info(this)

    class(t_particle), intent(in) :: this

    print "(5x,a)", "============================="
    print "(6x,a,1x,f3.0,1x,a)", 'Particle from:', this%originNum, this%originName
    print "(6x,a,1x,f10.6)", 'Lon:', this%xPos
    print "(6x,a,1x,f10.6)", 'Lat:', this%yPos
    print "(6x,a,1x,f10.6)", 'Depth [m]:', this%zPos
    print "(6x,a,1x,f10.6)", 'Age:', this%age
    print "(6x,a,1x,l1)", 'Active:', this%isActive
    if (.not. this%isActive) then
      print "(6x,a,1x,f3.0)", 'State:', this%state
    end if
    print "(5x,a)", "============================="

    return
  end subroutine
  !===========================================
  subroutine update(this)

    class(t_particle), intent(inout) :: this

    this%xPos = xnew
    this%yPos = ynew
    this%zPos = znew ! This will always stay the same if run_3d=.false.

    this%u = 0.5 * (pvelu + pvelunew)
    this%v = 0.5 * (pvelv + pvelvnew)
    this%w = 0.5 * (pvelw + pvelwnew)

    this%trajLen = this%trajLen + &
                   sqrt((this%u * dt)**2 + &
                        (this%v * dt)**2 + &
                        (this%w * dt)**2) ! This will always be 0 if run_3d=.false.

    this%age = this%age + dt

    return
  end subroutine update
  !===========================================
  subroutine check_age(this, max_age)

    class(t_particle), intent(inout)  :: this
    real(rk), intent(in)            :: max_age

    if (this%age > max_age) then
      this%isActive = .false.
    end if

  end subroutine check_age
  !===========================================
  subroutine check_depth(this, zaxarr, i, j, zval, change_state)

    class(t_particle), intent(inout) :: this
    logical, intent(in)            :: change_state
    integer, intent(in)            :: i, j
    real(rk), intent(in)           :: zaxarr(nx, ny, nlevels)
    real(rk)                       :: elev
    real(rk), intent(inout)        :: zval

    elev = sealevel(zaxarr, ig, jg)
    if (zval > elev) then
      DBG, "Setting particle to sealevel"
      zval = elev
    end if
    if (zval < -1.0 * depdata(i, j)) then
      DBG, "Setting particle to bottom depth"
      zval = -1.0 * depdata(i, j)
      if (change_state) then
        this%state = 3
        if (this%kill_bdy) this%isActive = .false.
      end if
    end if

  end subroutine check_depth
  !===========================================
  subroutine check_beached_bdy(this)

    class(t_particle), intent(inout) :: this

    dbghead(check_beached_bdy)

    call get_indices2d(this%xPos, this%yPos, x0, y0, dx, dy, ig, jg, igr, jgr)
    if (seamask(ig, jg) .eq. 3) then
      this%timeOnBeach = this%timeOnBeach + dt
    end if
    if (seamask(ig, jg) .eq. 2) then
      this%timeOnBeach = 0.d0
    end if
    !---------------------------------------------
    ! Change state if beaching time exceeded or on boundary
    if (this%timeOnBeach .ge. this%beachingtime) then
      if (this%kill_bch) this%isActive = .false.
      this%state = 1.
    end if

    if (seamask(ig, jg) .eq. 4) then
      if (this%kill_bdy) this%isActive = .false.
      this%state = 2.
    end if

    dbgtail(check_beached_bdy)
    return
  end subroutine check_beached_bdy

end module particle_type
!===================================================
module particle_vars
  !----------------------------------------------------------------
  ! This module includes variables related to particles:
  ! - number of particles, initial locations or something (maybe)...
  ! - anything else?
  ! TODO: Initial coordinates from netCDF
  !----------------------------------------------------------------
  use precdefs
  use errors
  use particle_type
  use nc_manager, only: nc_read_real_1d, nc_read_real_2d, nc_get_dim, nc_var_exists
  use domain_vars, only: lons, lats, nx, ny, seamask, x0, y0, dx, dy
  use fields, only: get_indices2d
  use time_vars, only: nTimes

  logical                     :: kill_beached, kill_boundary ! Set isActive=.false. if beached or on boundary?
  integer                     :: inputstep, &                ! How often are particles released?
                                 particle_init_method, &     ! Read initial positions (1 - txt, 2 - .nc)
                                 nParticles, &               ! Number of particles
                                 nInitTimes, &
                                 runparts = 0                ! Number of particles to loop over
  real(rk)                    :: max_age                     ! Lifetime (for all particles) in timesteps
  character(len=256)          :: coordfile                   ! File containing particle locations at init.
  type(t_particle), allocatable :: particles(:)                ! Array of particles
  ! real(rk), allocatable       :: initCoords(:, :)            ! Array of initial coordinates

  type, private :: t_initial_position
    integer                             :: next_idx = 1, &
                                           time_idx = 1, &
                                           n_init_particles
    real(rk), allocatable, dimension(:) :: x, y, z, &
                                           rho, radius, &
                                           beaching_time, &
                                           id
  contains
    procedure :: allocate_n_init_particles
    procedure :: check_initial_coordinates
  end type t_initial_position

  type(t_initial_position), allocatable :: initCoords(:)
!===================================================
contains
!===========================================
  subroutine allocate_n_init_particles(this)
    class(t_initial_position), intent(inout) :: this
    integer :: n_init_p

    n_init_p = this%n_init_particles

    allocate (this%x(n_init_p), this%y(n_init_p), this%z(n_init_p), &
              this%rho(n_init_p), this%radius(n_init_p), &
              this%beaching_time(n_init_p), this%id(n_init_p))

    this%x = 0.0d0; this%y = 0.0d0; this%z = 0.0d0; 
    this%rho = 0.0d0; this%radius = 0.0d0; 
    this%beaching_time = 0.0d0; this%id = 0.0d0

  end subroutine allocate_n_init_particles
  !===========================================
  subroutine init_particles_from_coordfile
    !---------------------------------------------
    ! Allocate array for estimated amount of particles
    !---------------------------------------------
    integer :: ipart

    !print
    FMT1, "======== Init particles ========"

    allocate (initCoords(1))

    open (COORDFILE, file=trim(coordfile), action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error("init_particles_from_coordfile", "Failed to open "//trim(coordfile), ierr)
    read (COORDFILE, *) initCoords(1)%n_init_particles
    call initCoords(1)%allocate_n_init_particles
    do ipart = 1, initCoords(1)%n_init_particles
      read (COORDFILE, *, iostat=ierr) initCoords(1)%x(ipart), initCoords(1)%y(ipart), &
        initCoords(1)%z(ipart), initCoords(1)%id(ipart), &
        initCoords(1)%beaching_time(ipart), initCoords(1)%rho(ipart), initCoords(1)%radius(ipart)
      if (ierr .ne. 0) call throw_error("init_particles_from_coordfile", "Failed to read from "//trim(coordfile), ierr)
    end do
    close (COORDFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error("init_particles_from_coordfile", "Failed to close "//trim(coordfile), ierr)

    call initCoords(1)%check_initial_coordinates

    nParticles = initCoords(1)%n_init_particles * (nTimes / inputstep) + initCoords(1)%n_init_particles
    allocate (particles(nParticles))
    FMT2, "Allocated array for", nParticles, "particles"

    !print
    FMT2, "Finished init particles"

  end subroutine init_particles_from_coordfile
  !===========================================
  subroutine init_particles_from_netcdf
    !---------------------------------------------
    ! Allocate array for estimated amount of particles
    !---------------------------------------------
    integer :: itime
    real(rk), allocatable :: nInitParticles(:)

    FMT1, "======== Init particles ========"

    call nc_get_dim(trim(coordfile), "time", nInitTimes)
    allocate (initCoords(nInitTimes), nInitParticles(nInitTimes))

    call nc_read_real_1d(trim(coordfile), "n_particles", nInitTimes, nInitParticles)

    debug(nInitTimes); debug(nInitParticles)

    do itime = 1, nInitTimes
      initCoords(itime)%time_idx = itime
      if (itime < nInitTimes) then
        initCoords(itime)%next_idx = itime + 1
      else
        ! Could also be periodic (last next_idx = 1)
        initCoords(itime)%next_idx = itime
      end if

      debug(itime)
      debug(int(nInitParticles(itime)))

      initCoords(itime)%n_init_particles = int(nInitParticles(itime))
      call initCoords(itime)%allocate_n_init_particles

      if (nc_var_exists(trim(coordfile), "x")) then
        call nc_read_real_2d(trim(coordfile), "x", 1, int(nInitParticles(itime)), initCoords(itime)%x, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        initCoords(itime)%x = 0.0d0
      end if

      if (nc_var_exists(trim(coordfile), "y")) then
        call nc_read_real_2d(trim(coordfile), "y", 1, int(nInitParticles(itime)), initCoords(itime)%y, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        initCoords(itime)%y = 0.0d0
      end if

      if (nc_var_exists(trim(coordfile), "z")) then
        call nc_read_real_2d(trim(coordfile), "z", 1, int(nInitParticles(itime)), initCoords(itime)%z, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        initCoords(itime)%z = 0.0d0
      end if

      if (nc_var_exists(trim(coordfile), "id")) then
        call nc_read_real_2d(trim(coordfile), "id", 1, int(nInitParticles(itime)), initCoords(itime)%id, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        initCoords(itime)%id = 0.0d0
      end if

      if (nc_var_exists(trim(coordfile), "beaching_time")) then
        call nc_read_real_2d(trim(coordfile), "beaching_time", 1, int(nInitParticles(itime)), initCoords(itime)%beaching_time, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        initCoords(itime)%beaching_time = 86400.0d0
      end if

      if (nc_var_exists(trim(coordfile), "rho")) then
        call nc_read_real_2d(trim(coordfile), "rho", 1, int(nInitParticles(itime)), initCoords(itime)%rho, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        initCoords(itime)%rho = 30.0d0
      end if

      if (nc_var_exists(trim(coordfile), "radius")) then
        call nc_read_real_2d(trim(coordfile), "radius", 1, int(nInitParticles(itime)), initCoords(itime)%radius, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        initCoords(itime)%radius = 0.001
      end if

      call initCoords(itime)%check_initial_coordinates

    end do

    if (nInitTimes > 1) then
      nParticles = int(sum(nInitParticles))
    else
      nParticles = int(nInitParticles(1)) * (nTimes / inputstep) + int(nInitParticles(1))
    end if
    allocate (particles(nParticles))
    FMT2, "Allocated array for", nParticles, "particles"

    FMT2, "Finished init particles"

  end subroutine init_particles_from_netcdf
  !===========================================
  subroutine check_initial_coordinates(this)
    class(t_initial_position), intent(in) :: this
    integer :: ipart, i, j, on_land = 0

    do ipart = 1, this%n_init_particles
      if (this%x(ipart) < lons(1) .or. this%x(ipart) > lons(nx) .or. &
          this%y(ipart) < lats(1) .or. this%y(ipart) > lats(ny)) then
        ERROR, "Particle", ipart, ":", &
          this%x(ipart), this%y(ipart)
        call throw_error("check_initial_coordinates", "Particle initialised outside of domain")
      end if
      call get_indices2d(this%x(ipart), this%y(ipart), x0, y0, dx, dy, i, j)
      if (seamask(i, j) == LAND) then
        ! call throw_warning("check_initial_coordinates", "Particle initialised on land")
        on_land = on_land + 1
      end if
    end do
    if (on_land > 0) then
      call throw_warning("check_initial_coordinates", "Particles initialised on land")
      WARNING, on_land, "particles on land, time idx = ", this%time_idx
      ! else
      !   FMT2, "Initial coordinates OK"
    end if

    return
  end subroutine check_initial_coordinates

end module particle_vars
