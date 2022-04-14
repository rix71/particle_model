#include "cppdefs.h"
module mod_physics
  !----------------------------------------------------------------
  ! This module contains the physics methods
  ! TODO:
  ! - diffusion
  ! - biofouling
  ! - drag
  !----------------------------------------------------------------
  use mod_errors
  use mod_interp
  use mod_precdefs
  use mod_params
  use time_vars, only: dt
  use mod_particle, only: t_particle
  use mod_fieldset, only: t_fieldset
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: vertical_velocity, init_diffusion, diffuse, update_Ah_Smagorinsky_full_field
  !---------------------------------------------
  real(rk), allocatable :: Ah_field(:, :, :)
  !---------------------------------------------
  integer :: ierr
  !===================================================
contains
  !===========================================
  subroutine vertical_velocity(p, fieldset, time, density_method, viscosity_method)
    !---------------------------------------------
    ! Calculate the vertical velocity due to buoyancy
    ! TODO: Which timestep should be used? (original or t + dt?)
    !---------------------------------------------
    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)            :: time
    integer, intent(in)             :: density_method
    logical, intent(in)             :: viscosity_method
    real(rk)                        :: i, j, k
    real(rk)                        :: kooi_vel
    real(rk)                        :: rho = ONE ! Fluid density
    real(rk)                        :: T, S
    real(rk)                        :: mu  ! Dynamic viscosity
    real(rk)                        :: kin_visc

    dbghead(vertical_velocity)

    debug(time)
    debug(p%depth1)
    debug(p%w1)

    i = p%ir0
    j = p%jr0
    k = p%kr0

    ! Density
    select case (density_method)
    case (DENSITY)
      rho = fieldset%get("RHO", time, i, j, k)
    case (TEMP_SALT)
      T = fieldset%get("TEMP", time, i, j, k)
      S = fieldset%get("SALT", time, i, j, k)
      rho = seawater_density_from_temp_and_salt(T, S, p%depth0)
    case (DEFAULT_DENSITY)
      rho = rho_default
    end select

    ! Viscosity
    select case (viscosity_method)
    case (.true.)
      mu = fieldset%get("VISC", time, i, j, k)
    case (.false.)
      mu = mu_default
    end select
    kin_visc = mu / rho

    kooi_vel = Kooi_vertical_velocity(p%rho, p%radius, rho, kin_visc)

    p%depth1 = p%depth1 + (kooi_vel * dt)
    p%w1 = p%w1 + kooi_vel

    debug(p%depth1)
    debug(p%w1)

    dbgtail(vertical_velocity)
    return
  end subroutine vertical_velocity
  !===========================================
  real(rk) function Kooi_vertical_velocity(rho_p, rad_p, rho_env, kin_visc) result(res)
    !---------------------------------------------
    ! Calculate vertical velocity
    ! Reference: Kooi 2017
    !---------------------------------------------
    real(rk), intent(in) :: rho_p, rad_p, rho_env
    real(rk), intent(in) :: kin_visc  ! Kinematic viscosity
    real(rk)             :: d_star    ! Dimensionless diameter
    real(rk)             :: w_star    ! Dimensionless settling velocity
    real(rk)             :: delta_rho ! Density difference

    dbghead(Kooi_vertical_velocity)

    debug(rho_p); debug(rad_p); 
    debug(rho_env)

    ! kin_visc = mu_env / rho_env ! NOT USING VISCOSITY???
    delta_rho = rho_p - rho_env

    debug(kin_visc); 
    debug(delta_rho)

    d_star = (delta_rho * g * (2.*rad_p)**3.) / (rho_env * (kin_visc**2.)) ! g negative?
    if (d_star < 0.05) then
      DBG, "d_star < 0.05"
      w_star = 1.74e-4 * (d_star**2)
    else if (d_star > 5.e9) then
      DBG, "d_star > 5e9"
      w_star = 1000.
    else
      DBG, "0.05 > d_star > 5e9"
      w_star = 10.**(-3.76715 + (1.92944 * log10(d_star)) - (0.09815 * log10(d_star)**2.) &
                     - (0.00575 * log10(d_star)**3.) + (0.00056 * log10(d_star)**4.))
    end if

    debug(d_star); debug(w_star)

    if (delta_rho > ZERO) then
      DBG, "delta_rho > 0"; debug(delta_rho / rho_env)
      res = -1.0 * ((delta_rho / rho_env) * g * w_star * kin_visc)**(1./3.) ! Getting NaNs with -1*g
    else
      DBG, "delta_rho < 0"; debug(delta_rho / rho_env)
      res = (-1.0 * (delta_rho / rho_env) * g * w_star * kin_visc)**(1./3.)
    end if

    debug(res)

    dbgtail(Kooi_vertical_velocity)
    return
  end function Kooi_vertical_velocity
  !===========================================
  subroutine init_diffusion(nx, ny, nz)
    integer, intent(in) :: nx, ny, nz

    dbghead(init_diffusion)

    FMT1, "======== Init diffusion ========"

    FMT2, "Allocating field for Ah (nx, ny, nz): (", nx, ", ", ny, ", ", nz, ")"
    allocate (Ah_field(nx, ny, nz), stat=ierr)
    if (ierr .ne. 0) call throw_error("physics :: init_diffusion", "Could not allocate!")

    ! Random seed!!!

    dbgtail(init_diffusion)
    return
  end subroutine init_diffusion
  !===========================================
  subroutine update_Ah_Smagorinsky_full_field(fieldset, time, nx, ny, nz)
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in) :: time
    integer, intent(in) :: nx, ny, nz
    real(rk), dimension(nx, ny, nz) :: dudx, dudy, dvdx, dvdy
    real(rk) :: dx, dy
    ! real(rk) :: Ah_s_full(nx, ny, nz)

    dbghead(update_Ah_Smagorinsky_full_field)

    dudx = fieldset%get_gradient("U", time, 1)
    dudy = fieldset%get_gradient("U", time, 2)
    dvdx = fieldset%get_gradient("V", time, 1)
    dvdy = fieldset%get_gradient("V", time, 2)

    dx = fieldset%domain%dx
    dy = fieldset%domain%dy

    Ah_field = Cm_smagorinsky * (dx * dy) * sqrt(dudx**2 + dvdy**2 + HALF * (dudy + dvdx)**2)

    dbgtail(update_Ah_Smagorinsky_full_field)
    return
  end subroutine update_Ah_Smagorinsky_full_field
  !===========================================
  function get_Ah_Smagorinsky_full_field(fieldset, x, y, z) result(Ah_s)
    type(t_fieldset), intent(in)   :: fieldset
    real(rk) AIM_INTENT            :: x, y
    real(rk), optional, intent(in) :: z
    integer                        :: i, j, k
    real(rk)                       :: f1
    real(rk)                       :: x1, x2, &
                                      y1, y2, &
                                      z1, z2, &
                                      c11, c12, c21, c22, &
                                      c111, c121, c211, c221, c112, c122, c212, c222
    real(rk) :: Ah_s
    integer :: inbr
    integer :: seamask(fieldset%nx, fieldset%ny)

    dbghead(get_Ah_Smagorinsky_full_field)

    debug(x); debug(y)

    i = floor(x); debug(i)
    x1 = float(floor(x)); debug(x1)
    x2 = float(floor(x) + 1); debug(x2)
    if (x >= real(fieldset%nx, rk)) then
#ifdef SNAP_TO_BOUNDS
      ! Right boundary
      ! can only be trusted when the SNAP_TO_BOUNDS flag is defined
      DBG, "Snapping x down"
      i = fieldset%nx - 1; debug(i)
      x = real(fieldset%nx, rk); debug(x)
      x1 = real(fieldset%nx, rk) - ONE; debug(x1)
      x2 = real(fieldset%nx, rk); debug(x2)
#else
      call throw_error("physics :: get_Ah_Smagorinsky_full_field", "Index (i) out of bounds!")
#endif
    end if
    if (x < ONE) then
#ifdef SNAP_TO_BOUNDS
      ! Left boundary
      ! can only be trusted when the SNAP_TO_BOUNDS flag is defined
      DBG, "Snapping x up"
      i = 1; debug(i)
      x = ONE; debug(x)
      x1 = ONE; debug(x1)
      x2 = 2.*ONE; debug(x2)
#else
      call throw_error("physics :: get_Ah_Smagorinsky_full_field", "Index (i) out of bounds!")
#endif
    end if

    j = floor(y); debug(j)
    y1 = float(floor(y)); debug(y1)
    y2 = float(floor(y) + 1); debug(y2)
    if (y >= real(fieldset%ny, rk)) then
#ifdef SNAP_TO_BOUNDS
      ! Top boundary
      DBG, "Snapping y down"
      j = fieldset%ny - 1; debug(j)
      y = real(fieldset%ny, rk); debug(y)
      y1 = real(fieldset%ny, rk) - ONE; debug(y1)
      y2 = real(fieldset%ny, rk); debug(y2)
#else
      call throw_error("physics :: get_Ah_Smagorinsky_full_field", "Index (j) out of bounds!")
#endif
    end if
    if (y < ONE) then
#ifdef SNAP_TO_BOUNDS
      ! Bottom boundary
      DBG, "Snapping y up"
      j = 1; debug(j)
      y = ONE; debug(y)
      y1 = ONE; debug(y1)
      y2 = 2.*ONE; debug(y2)
#else
      call throw_error("physics :: get_Ah_Smagorinsky_full_field", "Index (j) out of bounds!")
#endif
    end if

    seamask = fieldset%domain%get_seamask()

    if (present(z)) then
      debug(z)
      DBG, "Spatial intepolation 3D"

      k = floor(z); debug(k)
      z1 = float(floor(z)); debug(z1)
      z2 = float(floor(z) + 1); debug(z2)

      if (seamask(i, j) == LAND) then
        DBG, "Nearest neighbour"
        f1 = Ah_field(i, j, k)
        do inbr = 1, 8
          if (seamask(i + nbrs(1, inbr), j + nbrs(2, inbr)) .ne. LAND) then
            f1 = Ah_field(i + nbrs(1, inbr), j + nbrs(2, inbr), k)
            DBG, "Stopped neighbour search"
            DBG, "Data from: ", i + nbrs(1, inbr), j + nbrs(2, inbr), k
            debug(inbr)
            exit
          end if
        end do
        debug(f1)

      else if (seamask(i, j) == BEACH) then
        DBG, "Single cell value"
        f1 = Ah_field(i, j, k); debug(f1)

      else
        DBG, "Trilinear interpolation"

        if (k == fieldset%nz) then
          ! Use the lower layer if the particle is exactly at the surface
          k = k - 1
          z1 = z1 - ONE
          z2 = z2 - ONE
        end if

        c111 = Ah_field(i, j, k); debug(c111)
        c121 = Ah_field(i, j + 1, k); debug(c121)
        c211 = Ah_field(i + 1, j, k); debug(c211)
        c221 = Ah_field(i + 1, j + 1, k); debug(c221)
        c112 = Ah_field(i, j, k + 1); debug(c112)
        c122 = Ah_field(i, j + 1, k + 1); debug(c122)
        c212 = Ah_field(i + 1, j, k + 1); debug(c212)
        c222 = Ah_field(i + 1, j + 1, k + 1); debug(c222)

        call trilinearinterp(x1, x2, y1, y2, z1, z2, c111, c121, c211, c221, c112, c122, c212, c222, x, y, z, f1)
        debug(f1)
      end if
    else
      DBG, "Spatial intepolation 2D"

      if (seamask(i, j) == LAND) then
        DBG, "Nearest neighbour"
        f1 = Ah_field(i, j, 1)
        do inbr = 1, 8
          if (seamask(i + nbrs(1, inbr), j + nbrs(2, inbr)) .ne. LAND) then
            f1 = Ah_field(i + nbrs(1, inbr), j + nbrs(2, inbr), 1)
            DBG, "Stopped neighbour search"
            DBG, "Data from: ", i + nbrs(1, inbr), j + nbrs(2, inbr), 1
            debug(inbr)
            exit
          end if
        end do
        debug(f1)

      else if (seamask(i, j) == BEACH) then
        DBG, "Single cell value"
        f1 = Ah_field(i, j, 1); debug(f1)

      else
        DBG, "Bilinear interpolation"

        c11 = Ah_field(i, j, 1); debug(c11)
        c12 = Ah_field(i, j + 1, 1); debug(c12)
        c21 = Ah_field(i + 1, j, 1); debug(c21)
        c22 = Ah_field(i + 1, j + 1, 1); debug(c22)

        call bilinearinterp(x1, x1, x2, x2, y1, y2, c11, c12, c21, c22, x, y, f1)
      end if
    end if

    Ah_s = f1
    debug(Ah_s)

    dbgtail(get_Ah_Smagorinsky_full_field)
    return
  end function get_Ah_Smagorinsky_full_field
  !===========================================
  subroutine Ah_Smagorinsky(fieldset, time, i, j, k, Ah_s)
    type(t_fieldset), intent(in)  :: fieldset
    real(rk), intent(in)          :: time
    integer, intent(in)           :: i, j
    integer, intent(in), optional :: k
    real(rk), intent(out)         :: Ah_s
    real(rk)                      :: dx, dy
    real(rk)                      :: i0, j0, i1, j1, ih, jh, kh
    real(rk)                      :: uh0, vh0, uh1, vh1, u0h, v0h, u1h, v1h
    real(rk)                      :: dudx, dudy, dvdx, dvdy

    dbghead(Ah_Smagorinsky)

    debug(time); debug(i); debug(j)

    i0 = real(i, rk); debug(i0)
    j0 = real(j, rk); debug(j0)

    i1 = real(i, rk) + ONE; debug(i1)
    j1 = real(j, rk) + ONE; debug(j1)

    ih = real(i, rk) + HALF; debug(ih)
    jh = real(j, rk) + HALF; debug(jh)

    if (present(k)) then
      debug(k)
      kh = real(k, rk) + HALF; debug(kh)

      uh0 = fieldset%get("U", time, ih, j0, kh); debug(uh0)
      vh0 = fieldset%get("V", time, ih, j0, kh); debug(vh0)

      uh1 = fieldset%get("U", time, ih, j1, kh); debug(uh1)
      vh1 = fieldset%get("V", time, ih, j1, kh); debug(vh1)

      u0h = fieldset%get("U", time, i0, jh, kh); debug(u0h)
      v0h = fieldset%get("V", time, i0, jh, kh); debug(v0h)

      u1h = fieldset%get("U", time, i1, jh, kh); debug(u1h)
      v1h = fieldset%get("V", time, i1, jh, kh); debug(v1h)
    else
      uh0 = fieldset%get("U", time, ih, j0); debug(uh0)
      vh0 = fieldset%get("V", time, ih, j0); debug(vh0)

      uh1 = fieldset%get("U", time, ih, j1); debug(uh1)
      vh1 = fieldset%get("V", time, ih, j1); debug(vh1)

      u0h = fieldset%get("U", time, i0, jh); debug(u0h)
      v0h = fieldset%get("V", time, i0, jh); debug(v0h)

      u1h = fieldset%get("U", time, i1, jh); debug(u1h)
      v1h = fieldset%get("V", time, i1, jh); debug(v1h)
    end if

    dx = fieldset%domain%dx
    dy = fieldset%domain%dy

    dudx = (u1h - u0h) / dx; debug(dudx)
    dudy = (uh1 - uh0) / dy; debug(dudy)
    dvdx = (v1h - v0h) / dx; debug(dvdx)
    dvdy = (vh1 - vh0) / dy; debug(dvdy)

    ! Ah_s = Cm_smagorinsky * sqrt(HALF * (dudx - dvdy)**2 + HALF * (dudy + dvdx)**2)
    Ah_s = Cm_smagorinsky * (dx * dy) * sqrt(dudx**2 + dvdy**2 + HALF * (dudy + dvdx)**2)

    debug(Ah_s)

    dbgtail(Ah_Smagorinsky)
    return
  end subroutine Ah_Smagorinsky
  !===========================================
  subroutine diffuse_2D(p, fieldset, time)

    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)            :: time
    real(rk)                        :: Ah
    real(rk)                        :: i, j
    real(rk)                        :: x0, x1, &
                                       y0, y1

    dbghead(diffuse_2D)

    i = p%ir1
    j = p%jr1

    call fieldset%domain%lonlat2xy(p%lon1, p%lat1, x0, y0)

    ! call Ah_Smagorinsky(fieldset, time, i, j, Ah_s=Ah)
    Ah = get_Ah_Smagorinsky_full_field(fieldset, i, j)

    x1 = x0 + normal_random() * sqrt(2 * Ah * dt)
    y1 = y0 + normal_random() * sqrt(2 * Ah * dt)

    call fieldset%domain%xy2lonlat(x1, y1, p%lon1, p%lat1)
    call fieldset%search_indices(x=p%lon1, y=p%lat1, i=p%i1, j=p%j1, ir=p%ir1, jr=p%jr1)

    dbgtail(diffuse_2D)
    return
  end subroutine diffuse_2D
  !===========================================
  subroutine diffuse_3D(p, fieldset, time)

    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in) :: time
    real(rk)                        :: Ah, kv
    real(rk)                        :: i, j, k
    real(rk)                        :: x0, x1, &
                                       y0, y1, &
                                       z0, z1

    dbghead(diffuse_3D)

    i = p%ir1
    j = p%jr1
    k = p%kr1

    call fieldset%domain%lonlat2xy(p%lon1, p%lat1, x0, y0)
    z0 = p%depth1

    ! call Ah_Smagorinsky(fieldset, time, i, j, k=k, Ah_s=Ah)
    Ah = get_Ah_Smagorinsky_full_field(fieldset, i, j, k)
    kv = diffusion_vert_const

    x1 = x0 + normal_random() * sqrt(2 * Ah * dt)
    y1 = y0 + normal_random() * sqrt(2 * Ah * dt)
#ifdef DIFFUSE_VERTICAL
    z1 = z0 + normal_random() * sqrt(2 * kv * dt)
#else
    z1 = z0
#endif

    call fieldset%domain%xy2lonlat(x1, y1, p%lon1, p%lat1)
    call fieldset%search_indices(x=p%lon1, y=p%lat1, i=p%i1, j=p%j1, ir=p%ir1, jr=p%jr1)

    dbgtail(diffuse_3D)
    return
  end subroutine diffuse_3D
  !===========================================
  subroutine diffuse(p, fieldset, time, dif_3d)

    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in) :: time
    logical, intent(in)             :: dif_3d

    select case (dif_3d)
    case (.true.)
      call diffuse_3D(p, fieldset, time)
    case (.false.)
      call diffuse_2D(p, fieldset, time)
    end select

    return
  end subroutine diffuse
  !===========================================
  real(rk) function seawater_density_from_temp_and_salt(T, S, Z)
    !---------------------------------------------
    ! Density of seawater at one atmosphere pressure
    ! Reference:  N.P. Fofonoff and R.C. Millard Jr.,1983,
    !             Unesco technical papers in marine science no. 44.
    !---------------------------------------------

    real(rk), intent(in) :: T, S, Z
    real(rk) :: SAu, CTu, Zu, deltaS
    real(rk) :: R000, R100, R200, R300, R400, &
                R500, R600, R010, R110, R210, &
                R310, R410, R510, R020, R120, &
                R220, R320, R420, R030, R130, &
                R230, R330, R040, R140, R240, &
                R050, R150, R060, R001, R101, &
                R201, R301, R401, R011, R111, &
                R211, R311, R021, R121, R221, &
                R031, R131, R041, R002, R102, &
                R202, R012, R112, R022, R003, &
                R103, R013
    real(rk) :: ss, tt, rz3, rz2, rz1, rz0, zz

    dbghead(seawater_density_from_temp_and_salt)

    debug(T); debug(S)

    SAu = 40 * 35.16504 / 35; CTu = 40; Zu = 1e4

    deltaS = 32

    R000 = 8.0189615746e+02; R100 = 8.6672408165e+02; R200 = -1.7864682637e+03; 
    R300 = 2.0375295546e+03; R400 = -1.2849161071e+03; R500 = 4.3227585684e+02; 
    R600 = -6.0579916612e+01; R010 = 2.6010145068e+01; R110 = -6.5281885265e+01; 
    R210 = 8.1770425108e+01; R310 = -5.6888046321e+01; R410 = 1.7681814114e+01; 
    R510 = -1.9193502195e+00; R020 = -3.7074170417e+01; R120 = 6.1548258127e+01; 
    R220 = -6.0362551501e+01; R320 = 2.9130021253e+01; R420 = -5.4723692739e+00; 
    R030 = 2.1661789529e+01; R130 = -3.3449108469e+01; R230 = 1.9717078466e+01; 
    R330 = -3.1742946532e+00; R040 = -8.3627885467e+00; R140 = 1.1311538584e+01; 
    R240 = -5.3563304045e+00; R050 = 5.4048723791e-01; R150 = 4.8169980163e-01; 
    R060 = -1.9083568888e-01; R001 = 1.9681925209e+01; R101 = -4.2549998214e+01; 
    R201 = 5.0774768218e+01; R301 = -3.0938076334e+01; R401 = 6.6051753097e+00; 
    R011 = -1.3336301113e+01; R111 = -4.4870114575e+00; R211 = 5.0042598061e+00; 
    R311 = -6.5399043664e-01; R021 = 6.7080479603e+00; R121 = 3.5063081279e+00; 
    R221 = -1.8795372996e+00; R031 = -2.4649669534e+00; R131 = -5.5077101279e-01; 
    R041 = 5.5927935970e-01; R002 = 2.0660924175e+00; R102 = -4.9527603989e+00; 
    R202 = 2.5019633244e+00; R012 = 2.0564311499e+00; R112 = -2.1311365518e-01; 
    R022 = -1.2419983026e+00; R003 = -2.3342758797e-02; R103 = -1.8507636718e-02; 
    R013 = 3.7969820455e-01; 
    !

    ss = sqrt((S + deltaS) / SAu); 
    tt = T / CTu; 
    zz = -Z / Zu; 
    rz3 = R013 * tt + R103 * ss + R003

    rz2 = (R022 * tt + R112 * ss + R012) * tt + (R202 * ss + R102) * ss + R002

    rz1 = (((R041 * tt + R131 * ss + R031) * tt + (R221 * ss + R121) * ss + R021) * tt &
           + ((R311 * ss + R211) * ss + R111) * ss + R011) * tt &
          + (((R401 * ss + R301) * ss + R201) * ss + R101) * ss + R001

    rz0 = (((((R060 * tt + R150 * ss + R050) * tt + (R240 * ss + R140) * ss + R040) * tt &
           + ((R330 * ss + R230) * ss + R130) * ss + R030) * tt + (((R420 * ss + R320) * ss + R220) * ss + R120) * ss + R020) * tt &
           + ((((R510 * ss + R410) * ss + R310) * ss + R210) * ss + R110) * ss + R010) * tt &
          + (((((R600 * ss + R500) * ss + R400) * ss + R300) * ss + R200) * ss + R100) * ss + R000

    seawater_density_from_temp_and_salt = ((rz3 * zz + rz2) * zz + rz1) * zz + rz0

    dbgtail(seawater_density_from_temp_and_salt)
    return
  end function seawater_density_from_temp_and_salt
  !===========================================
  real(rk) function normal_random() result(r)
    !---------------------------------------------
    ! Get normally distributed random number from
    ! uniform distribution given by 'call random_number'
    !---------------------------------------------
    real(rk) :: r1, r2

    call random_number(r1); r1 = 1 - r1
    call random_number(r2); r2 = 1 - r2

    r = sqrt(-2 * log(r1)) * cos(2 * pi * r2)

    return
  end function normal_random
end module mod_physics
