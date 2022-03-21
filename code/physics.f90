#include "cppdefs.h"
module mod_physics
  !----------------------------------------------------------------
  ! This module contains the physics methods
  ! TODO:
  ! - diffusion
  ! - biofouling
  ! - drag
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_params
  use time_vars, only: dt
  use mod_particle, only: t_particle
  use mod_fieldset, only: t_fieldset
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: vertical_velocity
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
  subroutine diffuse(xin, yin, zin, dif_3d)

    real(rk), intent(inout) :: xin, yin, zin
    logical, intent(in)     :: dif_3d

    xin = xin + normal_random() * sqrt(2 * Ah * dt)
    yin = yin + normal_random() * sqrt(2 * Ah * dt)
    if (dif_3d) zin = zin + normal_random() * sqrt(2 * kv * dt)

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
