#include "cppdefs.h"
module particle_type
  !----------------------------------------------------------------
  ! This is the particle type definition
  !----------------------------------------------------------------
  use precdefs
  use domain_vars, only: x0, y0, dx, dy, seamask
  use loop_vars, only: ig, jg, igr, jgr, xnew, ynew, znew, &
                       pvelu, pvelunew, pvelv, pvelvnew, pvelw, pvelwnew
  use time_vars, only: dt
  use field_vars, only: uspeed, vspeed, wspeed, &
                        uspeednew, vspeednew, wspeednew, &
                        run_3d
  use fields, only: get_indices2d
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: particle
  !---------------------------------------------
  ! Particle type
  type particle
    logical           :: isActive = .true.     ! Skip particle in loop if isActive == .false.
    integer           :: warnings = 0
    real              :: state = 0.            ! 0 - active, 1 - beached, 2 - on boundary
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
    real(rk)          :: beachingtime        ! Different particles may essentialy have different beaching times
    real(rk)          :: originNum           ! Origin of particle, number

  contains
    procedure :: print_particle_info
    procedure :: update
    procedure :: check_beached_bdy
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
    this%u = 0.5 * (pvelv + pvelvnew)
    this%w = 0.5 * (pvelw + pvelwnew)

    this%trajLen = this%trajLen + &
                   sqrt((this%u * dt)**2 + &
                        (this%v * dt)**2 + &
                        (this%w * dt)**2) ! This will always be 0 if run_3d=.false.

    this%age = this%age + dt

    return
  end subroutine update
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
      this%isActive = .false.
      this%state = 1.
    end if
    ! TODO: Check boundary from seamask
    if (ig .lt. 10) then
      DBG, "Setting isActive=.false."
      this%isActive = .false.
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
  !----------------------------------------------------------------
  use precdefs
  use particle_type

  integer                     :: inputstep       ! How often are particles released?
  integer                     :: nParticles      ! Number of particles
  integer                     :: nInitParticles  ! Number of particles released every "inputstep"
  integer                     :: runparts = 0    ! Number of particles to loop over
  character(len=256)          :: coordfile       ! File containing particle locations at init.
  type(particle), allocatable :: particles(:)    ! Array of particles
  real(rk), allocatable       :: initCoords(:, :) ! Array of initial coordinates

end module particle_vars
