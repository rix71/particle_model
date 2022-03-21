#include "cppdefs.h"
module mod_particle
  !----------------------------------------------------------------
  ! This is the particle type definition
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  ! use mod_domain_vars, only: x0, y0, dx, dy, seamask, depdata, nx, ny
  ! Pass loop vars into functions rather than import?
  use time_vars, only: dt
  use mod_params, only: run_3d
  use mod_domain, only: t_domain
  use mod_fieldset, only: t_fieldset
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_particle
  !---------------------------------------------
  ! Particle type
  type t_particle
    logical  :: is_active = .true.    ! Skip particle in loop if is_active == .false.
    logical  :: kill_bch, kill_bdy    ! Set is_active=.false. if beached or on boundary?
    integer  :: warnings = 0
    integer  :: state = ACTIVE        ! 0 - active, 1 - beached, 2 - on boundary, 3 - bottom
    integer  :: i0, j0, k0            ! Particle position (grid cell indices, original)
    real(rk) :: ir0, jr0, kr0         ! Particle position (real indices, original)
    integer  :: i1, j1, k1            ! Particle position (grid cell indices, t + dt)
    real(rk) :: ir1, jr1, kr1         ! Particle position (real indices, t + dt)
    real(rk) :: lon0 = ZERO          ! Particle position (original)
    real(rk) :: lat0 = ZERO          ! Particle position (original)
    real(rk) :: depth0 = ZERO        ! Particle position (original)
    real(rk) :: lon1 = ZERO          ! Particle position (t + dt)
    real(rk) :: lat1 = ZERO          ! Particle position (t + dt)
    real(rk) :: depth1 = ZERO        ! Particle position (t + dt)
    real(rk) :: u0 = ZERO            ! Particle velocity (original)
    real(rk) :: v0 = ZERO            ! Particle velocity (original)
    real(rk) :: w0 = ZERO            ! Particle velocity (original)
    real(rk) :: u1 = ZERO            ! Particle velocity (t + dt)
    real(rk) :: v1 = ZERO            ! Particle velocity (t + dt)
    real(rk) :: w1 = ZERO            ! Particle velocity (t + dt)
    real(rk) :: rho = ZERO           ! Particle density
    real(rk) :: radius = ZERO        ! Particle radius
    real(rk) :: age = ZERO           ! Particle age
    real(rk) :: max_age = ZERO       ! Particle maximum age
    real(rk) :: traj_len = ZERO      ! Particle trajectory length
    real(rk) :: time_on_beach = ZERO ! Time spent in the beach area
    real(rk) :: beaching_time         ! Different particles may essentialy have different beaching times
    real(rk) :: id                    ! Origin of particle, number

  contains
    procedure :: update
    procedure :: check_boundaries
    procedure :: check_age
    procedure :: check_depth
    procedure :: bounce
    procedure :: print_info
  end type t_particle

  interface t_particle
    module procedure :: ctor_particle
  end interface t_particle

  !===================================================
contains
  !===========================================
  type(t_particle) function ctor_particle(lon, lat, depth, id, beaching_time, rho, radius, max_age, kill_bch, kill_bdy, fieldset, time) result(p)
    real(rk), intent(in)         :: lon, lat, depth
    real(rk), intent(in)         :: id
    real(rk), intent(in)         :: beaching_time
    real(rk), intent(in)         :: rho, radius
    real(rk), intent(in)         :: max_age
    logical, intent(in)          :: kill_bch, kill_bdy
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in)         :: time

    p%lon0 = lon
    p%lat0 = lat
    p%depth0 = depth
    p%id = id
    p%beaching_time = beaching_time
    p%rho = rho
    p%radius = radius
    p%max_age = max_age
    p%kill_bch = kill_bch
    p%kill_bdy = kill_bdy

    call fieldset%search_indices(time, lon, lat, depth, i=p%i0, j=p%j0, k=p%k0, ir=p%ir0, jr=p%jr0, kr=p%kr0)

  end function ctor_particle
  !===========================================
  subroutine update(this, fieldset, time)

    class(t_particle), intent(inout) :: this
    class(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)             :: time

    this%age = this%age + dt
    call this%check_age()

    ! Only update the age if the particle is beached or otherwise not active (but still alive)
    if (this%state .ne. ACTIVE) return

    call this%check_boundaries(fieldset, time)

    this%i0 = this%i1
    this%j0 = this%j1
    this%k0 = this%k1
    this%ir0 = this%ir1
    this%jr0 = this%jr1
    this%kr0 = this%kr1

    this%lon0 = this%lon1
    this%lat0 = this%lat1
    this%depth0 = this%depth1 ! This will always stay the same if run_3d=.false.

    this%u0 = this%u1
    this%v0 = this%v1
    this%w0 = this%w1

    this%traj_len = this%traj_len + &
                    sqrt((this%u0 * dt)**2 + &
                         (this%v0 * dt)**2 + &
                         (this%w0 * dt)**2) ! This will always be 0 if run_3d=.false.

    return
  end subroutine update
  !===========================================
  subroutine check_age(this)

    class(t_particle), intent(inout) :: this

    if ((this%max_age > 0) .and. (this%age > this%max_age)) then
      this%is_active = .false.
    end if

  end subroutine check_age
  !===========================================
  subroutine check_depth(this, fieldset, t, change_state)
    !---------------------------------------------
    ! TODO: Interpolation for bathymetry?
    !---------------------------------------------

    class(t_particle), intent(inout) :: this
    class(t_fieldset), intent(in)    :: fieldset
    logical, intent(in)              :: change_state
    real(rk), intent(in)             :: t
    real(rk)                         :: dep
    real(rk)                         :: elev
    real(rk)                         :: zval

    dbghead(check_depth)

    zval = this%depth1
    elev = fieldset%sealevel(t, this%ir1, this%jr1)

    debug(zval); debug(elev)

    if (zval > elev) then
      DBG, "Setting particle to sealevel"
      zval = elev
    end if

    dep = fieldset%domain%get_bathymetry(this%ir1, this%jr1)
    debug(dep)

    if (zval < -1.0 * dep) then
      DBG, "Setting particle to bottom depth"
      zval = -1.0 * dep
      if (change_state) then
        this%state = BOTTOM
        if (this%kill_bdy) this%is_active = .false.
      end if
    end if
    this%depth1 = zval

    dbgtail(check_depth)
  end subroutine check_depth
  !===========================================
  subroutine bounce(this, fieldset)
    class(t_particle), intent(inout) :: this
    class(t_fieldset), intent(in) :: fieldset
    integer :: i0, i1, j0, j1
    real(rk) :: ir0, ir1, jr0, jr1
    integer :: di, dj
    real(rk) :: dx0, dx1, dxp
    real(rk) :: dy0, dy1, dyp

    dbghead(bounce)

    i0 = this%i0; debug(i0)
    i1 = this%i1; debug(i1)
    j0 = this%j0; debug(j0)
    j1 = this%j1; debug(j1)

    ir0 = this%ir0; debug(ir0)
    ir1 = this%ir1; debug(ir1)
    jr0 = this%jr0; debug(jr0)
    jr1 = this%jr1; debug(jr1)

    di = i1 - i0
    dj = j1 - j0

#ifdef DEBUG
    if ((abs(di) > 1) .or. (abs(dj) > 1)) then
      DBG, "|di| > 1"
      debug(di)
      debug(dj)
      call throw_warning("particle :: bounce", "di or dj greater than 1!")
      ! return ! then what?
    end if
#endif

    if (di > ZERO) then
      ! has moved right
      DBG, "RIGHT"
      dx0 = floor(ir1) - ir0
      dx1 = ir1 - floor(ir1)
      dxp = dx0 - dx1
      ir1 = ir0 + (dxp / abs(di))
      ! i1 = i1 - 1
      i1 = int(ir1)
    else if (di < ZERO) then
      ! has moved left
      DBG, "LEFT"
      dx0 = ir0 - floor(ir0)
      dx1 = floor(ir0) - ir1
      dxp = dx0 - dx1
      ir1 = ir0 - (dxp / abs(di))
      ! i1 = i1 + 1
      i1 = int(ir1)
    end if

    if (dj > ZERO) then
      ! has moved up
      DBG, "UP"
      dy0 = floor(jr1) - jr0
      dy1 = jr1 - floor(jr1)
      dyp = dy0 - dy1
      jr1 = jr0 + (dyp / abs(dj))
      ! j1 = j1 - 1
      j1 = int(jr1)
    else if (dj < ZERO) then
      ! has moved down
      DBG, "DOWN"
      dy0 = jr0 - floor(jr0)
      dy1 = floor(jr0) - jr1
      dyp = dy0 - dy1
      jr1 = jr0 - (dyp / abs(dj))
      ! j1 = j1 + 1
      j1 = int(jr1)
    end if

    this%i1 = i1; debug(i1)
    this%j1 = j1; debug(j1)
    this%ir1 = ir1; debug(ir1)
    this%jr1 = jr1; debug(jr1)

    this%lon1 = fieldset%domain%get_lons(ir1); debug(this%lon1)
    this%lat1 = fieldset%domain%get_lats(jr1); debug(this%lat1)

    dbgtail(bounce)
    return
  end subroutine bounce
  !===========================================
  subroutine check_boundaries(this, fieldset, time)

    class(t_particle), intent(inout) :: this
    class(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)             :: time
    integer                          :: i, j
    integer                          :: seamask_val
    integer                          :: n_bounce
    integer, parameter               :: max_bounce = 100
#ifdef DEBUG
    real(rk)                         :: lon_t0 = 0., lat_t0 = 0.
    real(rk)                         :: lon_t1 = 0., lat_t1 = 0.
    real(rk)                         :: lon_t2 = 0., lat_t2 = 0.
#endif

    dbghead(check_boundaries)

    ! call fieldset%domain%get_indices_2d(this%lon1, this%lat1, i=i, j=j)
    i = this%i1
    j = this%j1

    debug(i); debug(j)

    seamask_val = fieldset%domain%get_seamask(i=i, j=j)

    debug(seamask_val)

    select case (seamask_val)
    case (BEACH)
      DBG, "Case BEACH"
      this%time_on_beach = this%time_on_beach + dt
      !---------------------------------------------
      ! Change state if beaching time exceeded or on boundary
      if (this%time_on_beach >= this%beaching_time) then
        if (this%kill_bch) this%is_active = .false.
        this%state = BEACHED
      end if
    case (LAND)
      DBG, "Case LAND"
      n_bounce = 0
      do while (seamask_val == LAND)
        call this%bounce(fieldset)
        i = this%i1
        j = this%j1
        seamask_val = fieldset%domain%get_seamask(i, j)
        n_bounce = n_bounce + 1
        if (n_bounce > max_bounce) exit
      end do
      debug(n_bounce)

#ifdef DEBUG
      lon_t2 = this%lon1
      lat_t2 = this%lat1
      if ((lon_t0 == lon_t2) .and. (lat_t0 == lat_t2)) then
        DBG, "CIRCLE"
        debug(fieldset%domain%get_seamask(this%i1, this%j1))
      end if
      lon_t0 = lon_t1
      lat_t0 = lat_t1
      lon_t1 = lon_t2
      lat_t1 = lat_t2
#endif

      !---------------------------------------------
      ! The bounce can happen only in the beach zone, so add to time on beach
      this%time_on_beach = this%time_on_beach + dt
      !---------------------------------------------
      ! Change state if beaching time exceeded or on boundary
      if (this%time_on_beach >= this%beaching_time) then
        if (this%kill_bch) this%is_active = .false.
        this%state = BEACHED
      end if
    case (SEA)
      DBG, "Case SEA"
      this%time_on_beach = ZERO
    case (BOUNDARY)
      DBG, "Case BOUNDARY"
      if (this%kill_bdy) this%is_active = .false.
      this%state = ON_BOUNDARY
    end select

    if (run_3d) call this%check_depth(fieldset, time, .true.)

    dbgtail(check_boundaries)
    return
  end subroutine check_boundaries
  !===========================================
  subroutine print_info(this)

    class(t_particle), intent(in) :: this

    FMT1, "Indices"
    FMT2, this%i0, this%j0, this%k0
    FMT2, this%i1, this%j1, this%k1

    FMT1, "Real indices"
    FMT2, this%ir0, this%jr0, this%kr0
    FMT2, this%ir1, this%jr1, this%kr1

    FMT1, "Position (lon, lat, depth)"
    FMT2, this%lon0, this%lat0, this%depth0
    FMT2, this%lon1, this%lat1, this%depth1

    FMT1, "Velocity"
    FMT2, this%u0, this%v0, this%w0
    FMT2, this%u1, this%v1, this%w1

    FMT1, "Other characteristics"
    FMT2, "Age: ", this%age
    FMT2, "Distance: ", this%traj_len
    FMT2, "Radius: ", this%radius
    FMT2, "Density: ", this%rho
    FMT2, "Time on beach: ", this%time_on_beach
    FMT2, "Beaching time: ", this%beaching_time

  end subroutine print_info

end module mod_particle
!===================================================
module mod_particle_vars
  !----------------------------------------------------------------
  ! This module includes variables related to particles:
  ! - number of particles, initial locations or something (maybe)...
  ! - anything else?
  ! TODO: Initial coordinates from netCDF
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  use mod_particle
  use nc_manager, only: nc_read_real_1d, nc_read_real_2d, nc_get_dim, nc_var_exists
  use mod_domain_vars, only: domain
  use time_vars, only: nTimes
  implicit none
  !===================================================
  !---------------------------------------------
  logical                       :: kill_beached, kill_boundary ! Set is_active=.false. if beached or on boundary?
  integer                       :: inputstep, &                ! How often are particles released?
                                   particle_init_method, &     ! Read initial positions (1 - txt, 2 - .nc)
                                   n_particles, &              ! Number of particles
                                   n_init_times, &
                                   runparts = 0                ! Number of particles to loop over
  real(rk)                      :: max_age                     ! Lifetime (for all particles) in timesteps
  character(len=LEN_CHAR_L)     :: coordfile                   ! File containing particle locations at init.
  type(t_particle), allocatable :: particles(:)                ! Array of particles
  !---------------------------------------------
  type, private :: t_initial_position
    integer                             :: next_idx = 1, &
                                           time_idx = 1, &
                                           n_particles
    real(rk), allocatable, dimension(:) :: x, y, z, &
                                           rho, radius, &
                                           beaching_time, &
                                           id
  contains
    procedure :: allocate_n_init_particles
    procedure :: check_initial_coordinates
  end type t_initial_position
  !---------------------------------------------
  type(t_initial_position), allocatable :: init_coords(:)
  !---------------------------------------------
  integer :: ierr
  !===================================================
contains
  !===========================================
  subroutine allocate_n_init_particles(this)
    class(t_initial_position), intent(inout) :: this
    integer :: n_init_p

    n_init_p = this%n_particles

    allocate (this%x(n_init_p), this%y(n_init_p), this%z(n_init_p), &
              this%rho(n_init_p), this%radius(n_init_p), &
              this%beaching_time(n_init_p), this%id(n_init_p))

    this%x = ZERO; this%y = ZERO; this%z = ZERO; 
    this%rho = ZERO; this%radius = ZERO; 
    this%beaching_time = ZERO; this%id = ZERO

  end subroutine allocate_n_init_particles
  !===========================================
  subroutine check_initial_coordinates(this)
    class(t_initial_position), intent(in) :: this
    integer :: ipart, i, j, on_land = 0

    do ipart = 1, this%n_particles
      if ((this%x(ipart) < domain%get_lons(1)) .or. (this%x(ipart) > domain%get_lons(domain%nx)) .or. &
          (this%y(ipart) < domain%get_lats(1)) .or. (this%y(ipart) > domain%get_lats(domain%ny))) then
        ERROR, "Particle", ipart, ":", &
          this%x(ipart), this%y(ipart)
        call throw_error("particle_vars :: check_initial_coordinates", "Particle initialised outside of domain")
      end if
      call domain%get_indices_2d(this%x(ipart), this%y(ipart), i=i, j=j)
      if (domain%get_seamask(i, j) == LAND) then
        ! call throw_warning("particle_vars :: check_initial_coordinates", "Particle initialised on land")
        on_land = on_land + 1
      end if
    end do
    if (on_land > 0) then
      call throw_warning("particle_vars :: check_initial_coordinates", "Particles initialised on land")
      WARNING, on_land, "particles on land, time idx = ", this%time_idx
      ! else
      !   FMT2, "Initial coordinates OK"
    end if

    return
  end subroutine check_initial_coordinates
  !===========================================
  subroutine init_particles_from_coordfile
    !---------------------------------------------
    ! Allocate array for estimated amount of particles
    !---------------------------------------------
    integer :: ipart

    !print
    FMT1, "======== Init particles ========"

    allocate (init_coords(1))

    open (COORDFILE, file=trim(coordfile), action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error("particle_vars :: init_particles_from_coordfile", "Failed to open "//trim(coordfile), ierr)
    read (COORDFILE, *) init_coords(1)%n_particles
    call init_coords(1)%allocate_n_init_particles
    do ipart = 1, init_coords(1)%n_particles
      read (COORDFILE, *, iostat=ierr) init_coords(1)%x(ipart), init_coords(1)%y(ipart), &
        init_coords(1)%z(ipart), init_coords(1)%id(ipart), &
        init_coords(1)%beaching_time(ipart), init_coords(1)%rho(ipart), init_coords(1)%radius(ipart)
  if (ierr .ne. 0) call throw_error("particle_vars :: init_particles_from_coordfile", "Failed to read from "//trim(coordfile), ierr)
    end do
    close (COORDFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error("particle_vars :: init_particles_from_coordfile", "Failed to close "//trim(coordfile), ierr)

    call init_coords(1)%check_initial_coordinates()

    n_particles = init_coords(1)%n_particles * (nTimes / inputstep) + init_coords(1)%n_particles
    allocate (particles(n_particles))
    FMT2, "Allocated array for", n_particles, "particles"

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

    call nc_get_dim(trim(coordfile), "time", n_init_times)
    allocate (init_coords(n_init_times), nInitParticles(n_init_times))

    call nc_read_real_1d(trim(coordfile), "n_particles", n_init_times, nInitParticles)

    debug(n_init_times); debug(nInitParticles)

    do itime = 1, n_init_times
      init_coords(itime)%time_idx = itime
      if (itime < n_init_times) then
        init_coords(itime)%next_idx = itime + 1
      else
        ! Could also be periodic (last next_idx = 1)
        init_coords(itime)%next_idx = itime
      end if

      debug(itime)
      debug(int(nInitParticles(itime)))

      init_coords(itime)%n_particles = int(nInitParticles(itime))
      call init_coords(itime)%allocate_n_init_particles

      if (nc_var_exists(trim(coordfile), "x")) then
        call nc_read_real_2d(trim(coordfile), "x", 1, int(nInitParticles(itime)), init_coords(itime)%x, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        init_coords(itime)%x = ZERO
      end if

      if (nc_var_exists(trim(coordfile), "y")) then
        call nc_read_real_2d(trim(coordfile), "y", 1, int(nInitParticles(itime)), init_coords(itime)%y, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        init_coords(itime)%y = ZERO
      end if

      if (nc_var_exists(trim(coordfile), "z")) then
        call nc_read_real_2d(trim(coordfile), "z", 1, int(nInitParticles(itime)), init_coords(itime)%z, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        init_coords(itime)%z = ZERO
      end if

      if (nc_var_exists(trim(coordfile), "id")) then
        call nc_read_real_2d(trim(coordfile), "id", 1, int(nInitParticles(itime)), init_coords(itime)%id, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        init_coords(itime)%id = ZERO
      end if

      if (nc_var_exists(trim(coordfile), "beaching_time")) then
        call nc_read_real_2d(trim(coordfile), "beaching_time", 1, int(nInitParticles(itime)), init_coords(itime)%beaching_time, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        init_coords(itime)%beaching_time = 86400.0d0
      end if

      if (nc_var_exists(trim(coordfile), "rho")) then
        call nc_read_real_2d(trim(coordfile), "rho", 1, int(nInitParticles(itime)), init_coords(itime)%rho, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        init_coords(itime)%rho = 30.0d0
      end if

      if (nc_var_exists(trim(coordfile), "radius")) then
        call nc_read_real_2d(trim(coordfile), "radius", 1, int(nInitParticles(itime)), init_coords(itime)%radius, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        init_coords(itime)%radius = 0.001
      end if

      call init_coords(itime)%check_initial_coordinates()

    end do

    if (n_init_times > 1) then
      n_particles = int(sum(nInitParticles))
    else
      n_particles = int(nInitParticles(1)) * (nTimes / inputstep) + int(nInitParticles(1))
    end if
    allocate (particles(n_particles))
    FMT2, "Allocated array for", n_particles, "particles"

    FMT2, "Finished init particles"

  end subroutine init_particles_from_netcdf

end module mod_particle_vars
