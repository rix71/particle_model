#include "cppdefs.h"
#include "particle.h"
#include "field.h"
#include "file.h"
module mod_particle
  !----------------------------------------------------------------
  ! This is the particle type definition
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  ! use mod_domain_vars, only: x0, y0, dx, dy, seamask, depdata, nx, ny
  ! Pass loop vars into functions rather than import?
  use time_vars, only: dt
  use mod_params, only: run_3d, pi
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
    !---------------------------------------------
    logical  :: is_active = .true.          ! Skip particle in loop if is_active == .false.
    logical  :: kill_beached, kill_boundary ! Set is_active=.false. if beached or on boundary?
    integer  :: warnings = 0
    integer  :: state = ST_SUSPENDED        ! States: active, beached, on boundary, bottom. Enumerated in cppdefs.h
    real(rk) :: id                          ! Origin of particle, number
    !---------------------------------------------
    ! Indices
    integer  :: i0, j0, k0                  ! Position (grid cell indices, original)
    real(rk) :: ir0, jr0, kr0               ! Position (real indices, original)
    integer  :: i1, j1, k1                  ! Position (grid cell indices, t + dt)
    real(rk) :: ir1, jr1, kr1               ! Position (real indices, t + dt)
    !---------------------------------------------
    ! Coordinates
    real(rk) :: lon0 = ZERO                 ! Position (original)
    real(rk) :: lat0 = ZERO                 ! Position (original)
    real(rk) :: depth0 = ZERO               ! Position (original)
    real(rk) :: lon1 = ZERO                 ! Position (t + dt)
    real(rk) :: lat1 = ZERO                 ! Position (t + dt)
    real(rk) :: depth1 = ZERO               ! Position (t + dt)
    !---------------------------------------------
    ! Velocity
    real(rk) :: u0 = ZERO                   ! Velocity (original)
    real(rk) :: v0 = ZERO                   ! Velocity (original)
    real(rk) :: w0 = ZERO                   ! Velocity (original)
    real(rk) :: u1 = ZERO                   ! Velocity (t + dt)
    real(rk) :: v1 = ZERO                   ! Velocity (t + dt)
    real(rk) :: w1 = ZERO                   ! Velocity (t + dt)
    real(rk) :: vel_vertical = ZERO         ! Settling velocity (Kooi)
    !---------------------------------------------
    real(rk) :: rho = ZERO                  ! Density
    real(rk) :: rho0 = ZERO                 ! Initial density
    real(rk) :: radius = ZERO               ! Radius
    real(rk) :: radius0 = ZERO              ! Initial radius
    real(rk) :: age = ZERO                  ! Age
    real(rk) :: max_age = ZERO              ! Maximum age
    real(rk) :: traj_len = ZERO             ! Particle trajectory length
    real(rk) :: time_on_beach = ZERO        ! Time spent in the beach area
    real(rk) :: beaching_time               ! Different particles may essentialy have different beaching times
    !---------------------------------------------
    ! Environment variables
    real(rk) :: delta_rho = ZERO            ! Density difference between particle and surrounding water
    real(rk) :: kin_visc = ZERO              ! Kinematic viscosity surrounding particle
    real(rk) :: u_star = ZERO               ! Friction velocity surrounding particle
    !---------------------------------------------
    ! Biofouling variables
    real(rk) :: growth_biofilm = ZERO       ! Attached algal growth
    real(rk) :: h_biofilm = ZERO            ! Thickness of biofilm

  contains
    private
    procedure, public :: update
    procedure         :: check_boundaries
    procedure         :: check_age
    procedure, public :: check_depth
    procedure         :: bounce
    procedure         :: redirect
    procedure, public :: print_info
    procedure, public :: volume
    procedure, public :: surface_area
  end type t_particle

  interface t_particle
    module procedure :: ctor_particle
    module procedure :: ctor_particle_restart
  end interface t_particle

  !===================================================
contains
  !===========================================
  type(t_particle) function ctor_particle(lon, lat, depth, &
                                          id, beaching_time, &
                                          rho, radius, max_age, &
                                          kill_beached, kill_boundary, &
                                          fieldset, time) result(p)
    real(rk), intent(in)         :: lon, lat, depth
    real(rk), intent(in)         :: id
    real(rk), intent(in)         :: beaching_time
    real(rk), intent(in)         :: rho, radius
    real(rk), intent(in)         :: max_age
    logical, intent(in)          :: kill_beached, kill_boundary
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in)         :: time

    p%lon0 = lon
    p%lat0 = lat
    p%depth0 = depth
    p%id = id
    p%beaching_time = beaching_time
    p%rho = rho
    p%rho0 = rho
    p%radius = radius
    p%radius0 = radius
    p%max_age = max_age
    p%kill_beached = kill_beached
    p%kill_boundary = kill_boundary

    call fieldset%search_indices(time, lon, lat, depth, i=p%i0, j=p%j0, k=p%k0, ir=p%ir0, jr=p%jr0, kr=p%kr0)

  end function ctor_particle
  !===========================================
  type(t_particle) function ctor_particle_restart(lon, lat, depth, &
                                                  i0, j0, k0, &
                                                  ir0, jr0, kr0, &
                                                  id, beaching_time, &
                                                  rho, rho0, &
                                                  radius, radius0, &
                                                  h_biofilm, &
                                                  age, max_age, kill_beached, kill_boundary, &
                                                  u0, v0, w0, vel_vertical, &
                                                  traj_len, time_on_beach, is_active, state) result(p)
    real(rk), intent(in)         :: lon, lat, depth
    real(rk), intent(in)         :: id
    real(rk), intent(in)         :: beaching_time
    real(rk), intent(in)         :: rho, rho0
    real(rk), intent(in)         :: radius, radius0
    real(rk), intent(in)         :: h_biofilm
    real(rk), intent(in)         :: age, max_age
    real(rk), intent(in)         :: ir0, jr0, kr0
    real(rk), intent(in)         :: u0, v0, w0, vel_vertical
    real(rk), intent(in)         :: traj_len, time_on_beach
    integer, intent(in)          :: i0, j0, k0
    integer, intent(in)          :: state
    logical, intent(in)          :: kill_beached, kill_boundary, is_active

    p%lon0 = lon
    p%lat0 = lat
    p%depth0 = depth

    p%i0 = i0
    p%j0 = j0
    p%k0 = k0
    p%ir0 = ir0
    p%jr0 = jr0
    p%kr0 = kr0

    p%u0 = u0
    p%v0 = v0
    p%w0 = w0
    p%vel_vertical = vel_vertical

    p%traj_len = traj_len

    p%id = id
    p%is_active = is_active
    p%state = state

    p%beaching_time = beaching_time
    p%time_on_beach = time_on_beach

    p%rho = rho
    p%rho0 = rho0

    p%radius = radius
    p%radius0 = radius0
    p%h_biofilm = h_biofilm

    p%age = age
    p%max_age = max_age

    p%kill_beached = kill_beached
    p%kill_boundary = kill_boundary

  end function ctor_particle_restart
  !===========================================
  elemental real(rk) function volume(this)
    class(t_particle), intent(in) :: this

    volume = 4./3.*pi * this%radius0**3

  end function volume
  !===========================================
  elemental real(rk) function surface_area(this)
    class(t_particle), intent(in) :: this

    surface_area = 4.*pi * this%radius0**2.

  end function surface_area
  !===========================================
  subroutine update(this, fieldset, time)

    class(t_particle), intent(inout) :: this
    class(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)             :: time
    real(rk)                         :: x0, y0
    real(rk)                         :: x1, y1

    this%age = this%age + dt
    call this%check_age()

    ! Only update the age if the particle is beached or otherwise not active (but still alive)
    if (this%state < ST_SUSPENDED) return

    call this%check_boundaries(fieldset, time)

    call fieldset%domain%lonlat2xy(this%lon0, this%lat0, x0, y0)
    call fieldset%domain%lonlat2xy(this%lon1, this%lat1, x1, y1)
    this%traj_len = this%traj_len + &
                    sqrt((x1 - x0)**2 + &
                         (y1 - y0)**2 + &
                         (this%depth1 - this%depth0)**2) ! This will always be 0 if run_3d=.false.

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

    return
  end subroutine update
  !===========================================
  subroutine check_age(this)

    class(t_particle), intent(inout) :: this

    if ((this%max_age > ZERO) .and. (this%age > this%max_age)) then
      this%is_active = .false.
    end if

  end subroutine check_age
  !===========================================
  subroutine check_depth(this, fieldset, t)
    !---------------------------------------------
    ! TODO: Interpolation for bathymetry?
    !---------------------------------------------

    class(t_particle), intent(inout) :: this
    class(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)             :: t
    real(rk)                         :: dep
    real(rk)                         :: elev

#ifdef PARTICLE_SNAP_SEALVL
    elev = fieldset%sealevel(t, this%ir1, this%jr1)
    if (this%depth1 >= elev) then
      this%depth1 = elev
      return
    end if
#endif

    dep = fieldset%domain%get_bathymetry(this%i1, this%j1)

    ! The particle is past the bottom
    if (this%depth1 <= -1.0 * dep) then
      this%depth1 = -1.0 * dep
      this%state = ST_BOTTOM
      this%w1 = ZERO
      return
    end if

    ! Reset to SUSPENDED if resuspended
    if ((this%depth1 > -1.0 * dep) .and. (this%state == ST_BOTTOM)) then
      this%state = ST_SUSPENDED
      return
    end if

  end subroutine check_depth
  !===========================================
  subroutine bounce(this, fieldset)
    class(t_particle), intent(inout) :: this
    class(t_fieldset), intent(in) :: fieldset
    integer :: i0, i1, i_orig, &
               j0, j1, j_orig
    real(rk) :: ir0, ir1, ir_orig, &
                jr0, jr1, jr_orig
    real(rk) :: lon_orig, lat_orig
    integer :: di, dj
    real(rk) :: dx0, dx1, dxp
    real(rk) :: dy0, dy1, dyp
    integer :: n_iter, seamask_val
    integer, parameter :: max_iter = 10

    ! We can initialise seamask_val like this, since this is the reason we're in this subroutine in the first place
    seamask_val = DOM_LAND

    lon_orig = this%lon0
    lat_orig = this%lat0

    i0 = this%i0; 
    i1 = this%i1; 
    j0 = this%j0; 
    j1 = this%j1; 
    ir0 = this%ir0; 
    ir1 = this%ir1; 
    jr0 = this%jr0; 
    jr1 = this%jr1; 
    i_orig = i0
    j_orig = j0
    ir_orig = ir0
    jr_orig = jr0

    n_iter = 1
    do while (seamask_val == DOM_LAND)
      di = i1 - i0
      dj = j1 - j0

#ifdef DEBUG
      if ((abs(di) > 1) .or. (abs(dj) > 1)) then

        call throw_warning("particle :: bounce", "di or dj greater than 1!")
        ! return ! then what?
      end if
#endif

      if (di > ZERO) then
        ! has moved right

        dx0 = floor(ir1) - ir0
        dx1 = ir1 - floor(ir1)
        dxp = dx0 - dx1
        ir1 = ir0 + (dxp / abs(di))
        ! i1 = i1 - 1
        i1 = int(ir1)
      else if (di < ZERO) then
        ! has moved left

        dx0 = ir0 - floor(ir0)
        dx1 = floor(ir0) - ir1
        dxp = dx0 - dx1
        ir1 = ir0 - (dxp / abs(di))
        ! i1 = i1 + 1
        i1 = int(ir1)
      end if

      if (dj > ZERO) then
        ! has moved up

        dy0 = floor(jr1) - jr0
        dy1 = jr1 - floor(jr1)
        dyp = dy0 - dy1
        jr1 = jr0 + (dyp / abs(dj))
        ! j1 = j1 - 1
        j1 = int(jr1)
      else if (dj < ZERO) then
        ! has moved down

        dy0 = jr0 - floor(jr0)
        dy1 = floor(jr0) - jr1
        dyp = dy0 - dy1
        jr1 = jr0 - (dyp / abs(dj))
        ! j1 = j1 + 1
        j1 = int(jr1)
      end if

      seamask_val = fieldset%domain%get_seamask(i1, j1)

      if (n_iter > max_iter) then

        this%i1 = i_orig
        this%j1 = j_orig
        this%ir1 = ir_orig
        this%jr1 = jr_orig
        this%lon1 = lon_orig
        this%lat1 = lat_orig

        return
      end if
    end do

    this%i1 = i1; 
    this%j1 = j1; 
    this%ir1 = ir1; 
    this%jr1 = jr1; 
    this%lon1 = fieldset%domain%get_lons(ir1); 
    this%lat1 = fieldset%domain%get_lats(jr1); 
    return
  end subroutine bounce
  !===========================================
  subroutine redirect(this, fieldset)
    class(t_particle), intent(inout) :: this
    class(t_fieldset), intent(in) :: fieldset
    integer :: i0, i1, i2, i_orig, &
               j0, j1, j2, j_orig
    real(rk) :: ir0, ir1, ir2, ir_orig, &
                jr0, jr1, jr2, jr_orig
    real(rk) :: lon_orig, lat_orig
    integer :: di, dj
    real(rk) :: dir, djr
    real(rk) :: dx0, dx1, dxp
    real(rk) :: dy0, dy1, dyp
    real(rk) :: x0, x1, y0, y1
    integer :: rot
    integer :: seamask_val
    integer :: n_iter
    integer, parameter :: max_iter = 10

    lon_orig = this%lon0
    lat_orig = this%lat0

    call fieldset%domain%lonlat2xy(this%lon0, this%lat0, x0, y0)

    i0 = this%i0; 
    i1 = this%i1; 
    j0 = this%j0; 
    j1 = this%j1; 
    i2 = i0
    j2 = j0
    i_orig = i0
    j_orig = j0

    ir0 = this%ir0; 
    ir1 = this%ir1; 
    jr0 = this%jr0; 
    jr1 = this%jr1; 
    ir2 = ir0
    jr2 = jr0
    ir_orig = ir0
    jr_orig = jr0

    di = i1 - i0; 
    dj = j1 - j0; 
    dir = ir1 - ir0; 
    djr = jr1 - jr0; 
    seamask_val = fieldset%domain%get_seamask(i1, j1)

    n_iter = 1
    do while (seamask_val == DOM_LAND)

      ! Choosing the order based on displacement
      if (abs(dir) > abs(djr)) then
        call right_left()

        i1 = int(ir1); 
        j1 = int(jr1); 
        di = i1 - i0; 
        dj = j1 - j0; 
        dir = ir1 - ir0; 
        djr = jr1 - jr0; 
        seamask_val = fieldset%domain%get_seamask(i1, j1); 
        if (seamask_val /= DOM_LAND) then
          exit
        end if

        call up_down()

        i0 = i_orig
        j0 = j_orig
        ir0 = ir_orig
        jr0 = jr_orig

        i1 = int(ir1); 
        j1 = int(jr1); 
        di = i1 - i0; 
        dj = j1 - j0; 
        dir = ir1 - ir0; 
        djr = jr1 - jr0; 
      else
        call up_down()

        i1 = int(ir1); 
        j1 = int(jr1); 
        di = i1 - i0; 
        dj = j1 - j0; 
        dir = ir1 - ir0; 
        djr = jr1 - jr0; 
        seamask_val = fieldset%domain%get_seamask(i1, j1); 
        if (seamask_val /= DOM_LAND) then
          exit
        end if

        call right_left()

        i0 = i_orig
        j0 = j_orig
        ir0 = ir_orig
        jr0 = jr_orig

        i1 = int(ir1); 
        j1 = int(jr1); 
        di = i1 - i0; 
        dj = j1 - j0; 
        dir = ir1 - ir0; 
        djr = jr1 - jr0; 
      end if

      seamask_val = fieldset%domain%get_seamask(i1, j1); 
      n_iter = n_iter + 1
      if (n_iter > max_iter) then

        this%i1 = i_orig
        this%j1 = j_orig
        this%ir1 = ir_orig
        this%jr1 = jr_orig
        this%lon1 = lon_orig
        this%lat1 = lat_orig

        return
      end if
    end do

    x1 = x0 + fieldset%domain%dx * dir; 
    y1 = y0 + fieldset%domain%dy * djr; 
    call fieldset%domain%xy2lonlat(x1, y1, this%lon1, this%lat1)

    this%i1 = i1; 
    this%j1 = j1; 
    this%ir1 = ir1; 
    this%jr1 = jr1; 
    ! this%lon1 = fieldset%domain%get_lons(ir1);

    ! this%lat1 = fieldset%domain%get_lats(jr1);

    return

  contains
    !===========================================
    subroutine up_down()
      if (dj > ZERO) then
        ! Up

        if (dir > ZERO) then
          rot = 1
        else
          rot = -1
        end if

        dy0 = floor(jr1) - jr0; 
        dy1 = jr1 - floor(jr1); 
        dxp = abs((dir / djr) * dy1); 
        dxp = max(dxp, SMALL); 
        dyp = dy0 - min(dxp, dy0); 
        ir2 = ir0 + dir + (dxp * rot); 
        jr2 = jr0 + dyp; 
        i0 = i1
        j0 = j1
        ir0 = ir1
        jr0 = jr1
        ir1 = ir2
        jr1 = jr2
      else if (dj < ZERO) then
        ! Down

        if (dir < ZERO) then
          rot = -1
        else
          rot = 1
        end if

        dy0 = jr0 - floor(jr0); 
        dy1 = floor(jr0) - jr1; 
        dxp = abs((dir / djr) * dy1); 
        dxp = max(dxp, SMALL); 
        dyp = dy0 - min(dxp, dy0); 
        ir2 = ir0 + dir + (dxp * rot); 
        jr2 = jr0 - dyp; 
        i0 = i1
        j0 = j1
        ir0 = ir1
        jr0 = jr1
        ir1 = ir2
        jr1 = jr2
      end if
    end subroutine up_down
    !===========================================
    subroutine right_left()
      if (di > ZERO) then
        ! Right

        if (djr < ZERO) then
          rot = -1
        else
          rot = 1
        end if

        dx0 = floor(ir1) - ir0; 
        dx1 = ir1 - floor(ir1); 
        dyp = abs((djr / dir) * dx1); 
        dyp = max(dyp, SMALL); 
        dxp = dx0 - min(dyp, dx0); 
        ir2 = ir0 + dxp; 
        jr2 = jr0 + djr + (dyp * rot); 
        i0 = i1
        j0 = j1
        ir0 = ir1
        jr0 = jr1
        ir1 = ir2
        jr1 = jr2
      else if (di < ZERO) then
        ! Left

        if (djr > ZERO) then
          rot = 1
        else
          rot = -1
        end if

        dx0 = ir0 - floor(ir0); 
        dx1 = floor(ir0) - ir1; 
        dyp = abs((djr / dir) * dx1); 
        dyp = max(dyp, SMALL); 
        dxp = dx0 - min(dyp, dx0); 
        ir2 = ir0 - dxp; 
        jr2 = jr0 + djr + (dyp * rot); 
        i0 = i1
        j0 = j1
        ir0 = ir1
        jr0 = jr1
        ir1 = ir2
        jr1 = jr2
      end if
    end subroutine right_left

  end subroutine redirect
  !===========================================
  subroutine check_boundaries(this, fieldset, time)

    class(t_particle), intent(inout) :: this
    class(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)             :: time
    integer                          :: i, j
    integer                          :: seamask_val
#ifdef DEBUG
    real(rk)                         :: lon_t0 = 0., lat_t0 = 0.
    real(rk)                         :: lon_t1 = 0., lat_t1 = 0.
    real(rk)                         :: lon_t2 = 0., lat_t2 = 0.
#endif

    ! call fieldset%domain%get_indices_2d(this%lon1, this%lat1, i=i, j=j)
    i = this%i1
    j = this%j1

    seamask_val = fieldset%domain%get_seamask(i=i, j=j)

    select case (seamask_val)
    case (DOM_BEACH)

      this%time_on_beach = this%time_on_beach + dt
      !---------------------------------------------
      ! Change state if beaching time exceeded or on boundary
#ifndef PARTICLE_BEACH_IMMEDIATELY
      if (this%time_on_beach >= this%beaching_time) then
#endif
        if (this%kill_beached) this%is_active = .false.
        this%state = ST_BEACHED
#ifndef PARTICLE_BEACH_IMMEDIATELY
      end if
#endif
    case (DOM_LAND)

#if defined PARTICLE_BOUNCE
      call this%bounce(fieldset)
#elif defined PARTICLE_REDIRECT
      call this%redirect(fieldset)
#elif defined PARTICLE_BEACH_IMMEDIATELY
      if (this%kill_beached) this%is_active = .false.
      this%state = ST_BEACHED
#else
      this%i1 = this%i0
      this%j1 = this%j0
      this%ir1 = this%ir0
      this%jr1 = this%jr0
      this%lon1 = this%lon0
      this%lat1 = this%lat0
      this%depth1 = this%depth0
#endif

#if defined DEBUG && (defined PARTICLE_BOUNCE || defined PARTICLE_REDIRECT)
      lon_t2 = this%lon1
      lat_t2 = this%lat1
      if ((lon_t0 == lon_t2) .and. (lat_t0 == lat_t2)) then

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
        if (this%kill_beached) this%is_active = .false.
        this%state = ST_BEACHED
      end if
    case (DOM_SEA)

      this%time_on_beach = ZERO
    case (DOM_BOUNDARY)

      if (this%kill_boundary) this%is_active = .false.
      this%state = ST_BOUNDARY
    end select

    if (run_3d) call this%check_depth(fieldset, time)

    return
  end subroutine check_boundaries
  !===========================================
  subroutine print_info(this)

    class(t_particle), intent(in) :: this

    FMT1, "Indices (i, j, k)"
    FMT2, this%i0, this%j0, this%k0
    FMT2, this%i1, this%j1, this%k1

    FMT1, "Real indices (ir, jr, kr)"
    FMT2, this%ir0, this%jr0, this%kr0
    FMT2, this%ir1, this%jr1, this%kr1

    FMT1, "Position (lon, lat, depth)"
    FMT2, this%lon0, this%lat0, this%depth0
    FMT2, this%lon1, this%lat1, this%depth1

    FMT1, "Velocity (u, v, w)"
    FMT2, this%u0, this%v0, this%w0
    FMT2, this%u1, this%v1, this%w1

    FMT1, "Velocity from buoyancy"
    FMT2, this%vel_vertical

    FMT1, "State"
    FMT2, "state: ", this%state
    FMT2, "is_active: ", this%is_active

    FMT1, "Other characteristics"
    FMT2, "Age: ", this%age
    FMT2, "Distance: ", this%traj_len
    FMT2, "Radius (R, R0): ", this%radius, this%radius0
    FMT2, "Density (rho, rho0): ", this%rho, this%rho0
    FMT2, "Time on beach: ", this%time_on_beach
    FMT2, "Beaching time: ", this%beaching_time
    FMT2, "Biofilm thickness: ", this%h_biofilm

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
  use time_vars, only: nTimes, run_start_dt, dt
  use run_params, only: runid, restart, restart_path
  use mod_datetime, only: t_datetime, datetime_from_netcdf
  use mod_fieldset, only: t_fieldset
  implicit none
  !===================================================
  !---------------------------------------------
  logical                       :: kill_beached, kill_boundary ! Set is_active=.false. if beached or on boundary?
  integer                       :: inputstep, &                ! How often are particles released?
                                   particle_init_method, &     ! Read initial positions (1 - txt, 2 - .nc)
                                   n_particles, &              ! Number of particles
                                   n_restart_particles, &
                                   n_init_times, &
                                   runparts = 0, &             ! Number of particles to loop over
                                   i_release = 1
  real(rk)                      :: max_age                     ! Lifetime (for all particles) in timesteps
  character(len=LEN_CHAR_L)     :: coordfile                   ! File containing particle locations at init.
  type(t_particle), allocatable :: particles(:)                ! Array of particles
  !---------------------------------------------
  type, private :: t_initial_position
    type(t_datetime)                    :: release_date
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
  subroutine init_particles()

    select case (particle_init_method)
    case (TXT_FILE)
      call init_particles_from_coordfile ! particle.f90
    case (NC_FILE)
      call init_particles_from_netcdf    ! particle.f90
    end select

    if (restart) then
      call check_restart_file
      allocate (particles(n_particles + n_restart_particles))
      FMT2, "Allocated array for", n_particles + n_restart_particles, "particles"
      call read_restart_file
    else
      allocate (particles(n_particles))
      FMT2, "Allocated array for", n_particles, "particles"
    end if

  end subroutine init_particles
  !===========================================
  subroutine check_restart_file()

    character(len=LEN_CHAR_L) :: restart_filename
    character(len=14) :: time_str
    logical  :: file_exists

    write (time_str, '(i0.14)') run_start_dt%shortDate(include_time=.true.)

    restart_filename = trim(restart_path)//"/"//trim(runid)//"."//trim(time_str)//".restart.dat"

    inquire (file=trim(restart_filename), exist=file_exists)
    if (file_exists) then
      open (RESTARTFILE, file=trim(restart_filename), action='read', iostat=ierr)
      read (RESTARTFILE, *) n_restart_particles
      close (RESTARTFILE)
    else
      call throw_error("particle :: read_restart_file", "No restart file found.")
    end if

  end subroutine check_restart_file
  !===========================================
  subroutine read_restart_file()
    !---------------------------------------------
    ! Read restart file (latest position)
    ! Initialise in particles array
    ! Update runparts
    !---------------------------------------------
    character(len=LEN_CHAR_L) :: restart_filename
    character(len=14) :: time_str
    integer  :: i
    real(rk) :: lon, lat, depth
    real(rk) :: id
    real(rk) :: beaching_time
    real(rk) :: rho, rho0
    real(rk) :: radius, radius0
    real(rk) :: h_biofilm
    real(rk) :: age, max_age
    real(rk) :: ir0, jr0, kr0
    real(rk) :: u0, v0, w0, vel_vertical
    real(rk) :: traj_len, time_on_beach
    integer  :: i0, j0, k0
    integer  :: state
    logical  :: kill_beached, kill_boundary, is_active

    write (time_str, '(i0.14)') run_start_dt%shortDate(include_time=.true.)

    restart_filename = trim(restart_path)//"/"//trim(runid)//"."//trim(time_str)//".restart.dat"

    open (RESTARTFILE, file=trim(restart_filename), action='read', status='old', iostat=ierr)
    read (RESTARTFILE, *)
    do i = 1, n_restart_particles
      read (RESTARTFILE, *) lon, lat, depth, &
        i0, j0, k0, &
        ir0, jr0, kr0, &
        id, beaching_time, &
        rho, rho0, &
        radius, radius0, &
        h_biofilm, &
        age, max_age, kill_beached, kill_boundary, &
        u0, v0, w0, vel_vertical, &
        traj_len, time_on_beach, is_active, state

      particles(i) = t_particle(lon, lat, depth, &
                                i0, j0, k0, &
                                ir0, jr0, kr0, &
                                id, beaching_time, &
                                rho, rho0, &
                                radius, radius0, &
                                h_biofilm, &
                                age, max_age, kill_beached, kill_boundary, &
                                u0, v0, w0, vel_vertical, &
                                traj_len, time_on_beach, is_active, state)
    end do
    close (RESTARTFILE)
    runparts = n_restart_particles

    return
  end subroutine read_restart_file
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
      if (domain%get_seamask(i, j) == DOM_LAND) then
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

    ! If inputstep is < 0, it means no additional particles will be released (for restart)
    if (inputstep > 0) then
      n_particles = init_coords(1)%n_particles * (nTimes / inputstep) + init_coords(1)%n_particles
    else
      n_particles = 0
    end if

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

    do itime = 1, n_init_times
      init_coords(itime)%time_idx = itime
      if (itime < n_init_times) then
        init_coords(itime)%next_idx = itime + 1
      else
        ! Could also be periodic (last next_idx = 1)
        ! will stop releasing particles when the init file runs out
        init_coords(itime)%next_idx = 1
        ! init_coords(itime)%next_idx = itime
      end if

      init_coords(itime)%n_particles = int(nInitParticles(itime))
      call init_coords(itime)%allocate_n_init_particles

      init_coords(itime)%release_date = datetime_from_netcdf(trim(coordfile), itime)

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

    do itime = 1, n_init_times
      if (run_start_dt <= init_coords(itime)%release_date) then
        i_release = itime
        call init_coords(itime)%release_date%print_short_date()
        exit
      end if
    end do

    if (inputstep > 0) then
      if (n_init_times > 1) then
        n_particles = int(sum(nInitParticles))
      else
        n_particles = int(nInitParticles(1)) * (nTimes / inputstep) + int(nInitParticles(1))
      end if
    else
      n_particles = 0
    end if

    FMT2, "Finished init particles"

  end subroutine init_particles_from_netcdf
  !===========================================
  subroutine release_particles(itime, date, fieldset, fieldset_time)
    integer, intent(in) :: itime
    type(t_datetime), intent(in) :: date
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in) :: fieldset_time
    integer :: ipart

    if (inputstep <= 0) return

    select case (particle_init_method)
    case (TXT_FILE)
      if (mod(itime, inputstep) /= 0) then
        return
      end if
    case (NC_FILE)
      if (date < init_coords(i_release)%release_date) then
        return
      end if
    end select

    FMT2, "Releasing ", init_coords(i_release)%n_particles, " new particles at itime = ", itime
    do ipart = 1, init_coords(i_release)%n_particles
      particles(ipart + runparts) = t_particle(lon=init_coords(i_release)%x(ipart), &
                                               lat=init_coords(i_release)%y(ipart), &
                                               depth=init_coords(i_release)%z(ipart), &
                                               id=init_coords(i_release)%id(ipart), &
                                               beaching_time=init_coords(i_release)%beaching_time(ipart), &
                                               rho=init_coords(i_release)%rho(ipart), &
                                               radius=init_coords(i_release)%radius(ipart), &
                                               max_age=max_age, &
                                               kill_beached=kill_beached, &
                                               kill_boundary=kill_boundary, &
                                               fieldset=fieldset, &
                                               time=fieldset_time)
    end do
    runparts = runparts + init_coords(i_release)%n_particles
    FMT2, runparts, "particles"
    i_release = init_coords(i_release)%next_idx

    return
  end subroutine release_particles

end module mod_particle_vars
