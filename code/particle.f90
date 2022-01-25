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
  public :: particle
  !---------------------------------------------
  ! Particle type
  type particle
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
  end type particle

  !===================================================
contains
  !===========================================
  subroutine print_particle_info(this)

    class(particle), intent(in) :: this

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

    class(particle), intent(inout) :: this

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

    class(particle), intent(inout)  :: this
    real(rk), intent(in)            :: max_age

    if (this%age > max_age) then
      this%isActive = .false.
    end if

  end subroutine check_age
  !===========================================
  subroutine check_depth(this, zaxarr, i, j, zval, change_state)

    class(particle), intent(inout) :: this
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

    class(particle), intent(inout) :: this

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
  use nc_manager, only: nc_read1d, nc_get_dim, nc_var_exists
  use domain_vars, only: lons, lats, nx, ny
  use time_vars, only: nTimes

  logical                     :: kill_beached, kill_boundary ! Set isActive=.false. if beached or on boundary?
  integer                     :: inputstep, &                ! How often are particles released?
                                 particle_init_method, &     ! Read initial positions (1 - txt, 2 - .nc)
                                 nParticles, &               ! Number of particles
                                 nInitParticles, &           ! Number of particles released every "inputstep"
                                 runparts = 0                ! Number of particles to loop over
  real(rk)                    :: max_age                     ! Lifetime (for all particles) in timesteps
  character(len=256)          :: coordfile                   ! File containing particle locations at init.
  type(particle), allocatable :: particles(:)                ! Array of particles
  ! real(rk), allocatable       :: initCoords(:, :)            ! Array of initial coordinates

  type, private :: initial_position
    real(rk), allocatable, dimension(:) :: x, y, z, &
                                           rho, radius, &
                                           beaching_time, &
                                           id
  contains
    procedure :: allocate_nInitParticles
  end type initial_position

  type(initial_position) :: initCoords
!===================================================
contains
!===========================================
  subroutine allocate_nInitParticles(this, n_init)
    class(initial_position), intent(inout) :: this
    integer, intent(in) :: n_init

    allocate (this%x(n_init), this%y(n_init), this%z(n_init), &
              this%rho(n_init), this%radius(n_init), &
              this%beaching_time(n_init), this%id(n_init))

    this%x = 0.0d0; this%y = 0.0d0; this%z = 0.0d0; 
    this%rho = 0.0d0; this%radius = 0.0d0; 
    this%beaching_time = 0.0d0; this%id = 0.0d0

  end subroutine allocate_nInitParticles
  !===========================================
  subroutine init_particles_from_coordfile
    !---------------------------------------------
    ! Allocate array for estimated amount of particles
    !---------------------------------------------
    integer :: ipart

    !print
    FMT1, "======== Init particles ========"

    open (COORDFILE, file=trim(coordfile), action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error("init_particles_from_coordfile", "Failed to open "//trim(coordfile), ierr)
    read (COORDFILE, *) nInitParticles
    ! allocate (initCoords(nInitParticles, 7))
    call initCoords%allocate_nInitParticles(nInitParticles)
    do ipart = 1, nInitParticles
      read (COORDFILE, *, iostat=ierr) initCoords%x(ipart), initCoords%y(ipart), &
        initCoords%z(ipart), initCoords%id(ipart), &
        initCoords%beaching_time(ipart), initCoords%rho(ipart), initCoords%radius(ipart)
      if (ierr .ne. 0) call throw_error("init_particles_from_coordfile", "Failed to read from "//trim(coordfile), ierr)
    end do
    close (COORDFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error("init_particles_from_coordfile", "Failed to close "//trim(coordfile), ierr)

    call check_initial_coordinates

    nParticles = nInitParticles * (nTimes / inputstep) + nInitParticles
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

    FMT1, "======== Init particles ========"

    call nc_get_dim(trim(coordfile), "particle", nInitParticles)

    call initCoords%allocate_nInitParticles(nInitParticles)

    if (nc_var_exists(trim(coordfile), "x")) then
      call nc_read1d(trim(coordfile), "x", nInitParticles, initCoords%x)
    else
      initCoords%x = 0.0d0
    end if

    if (nc_var_exists(trim(coordfile), "y")) then
      call nc_read1d(trim(coordfile), "y", nInitParticles, initCoords%y)
    else
      initCoords%y = 0.0d0
    end if

    if (nc_var_exists(trim(coordfile), "z")) then
      call nc_read1d(trim(coordfile), "z", nInitParticles, initCoords%z)
    else
      initCoords%z = 0.0d0
    end if

    if (nc_var_exists(trim(coordfile), "id")) then
      call nc_read1d(trim(coordfile), "id", nInitParticles, initCoords%id)
    else
      initCoords%id = 0.0d0
    end if

    if (nc_var_exists(trim(coordfile), "beaching_time")) then
      call nc_read1d(trim(coordfile), "beaching_time", nInitParticles, initCoords%beaching_time)
    else
      initCoords%beaching_time = 86400.0d0
    end if

    if (nc_var_exists(trim(coordfile), "rho")) then
      call nc_read1d(trim(coordfile), "rho", nInitParticles, initCoords%rho)
    else
      initCoords%rho = 30.0d0
    end if

    if (nc_var_exists(trim(coordfile), "radius")) then
      call nc_read1d(trim(coordfile), "radius", nInitParticles, initCoords%radius)
    else
      initCoords%radius = 0.001
    end if

    call check_initial_coordinates

    nParticles = nInitParticles * (nTimes / inputstep) + nInitParticles
    allocate (particles(nParticles))
    FMT2, "Allocated array for", nParticles, "particles"

    FMT2, "Finished init particles"

  end subroutine init_particles_from_netcdf
  !===========================================
  subroutine check_initial_coordinates

    integer :: ipart

    do ipart = 1, nInitParticles
      if (initCoords%x(ipart) < lons(1) .or. initCoords%x(ipart) > lons(nx) .or. &
          initCoords%y(ipart) < lats(1) .or. initCoords%y(ipart) > lats(ny)) then
        ERROR, "Particle", ipart, ":", &
          initCoords%x(ipart), initCoords%y(ipart)
        call throw_error("check_initial_coordinates", "Particle initialized outside of domain")
      end if
    end do
    FMT2, "Initial coordinates OK"

    return
  end subroutine check_initial_coordinates

end module particle_vars
