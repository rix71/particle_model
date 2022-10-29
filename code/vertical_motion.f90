#include "cppdefs.h"
module mod_vertical_motion
  !----------------------------------------------------------------
  ! Calculates the particles' vertical velocity
  !----------------------------------------------------------------
  use mod_errors
  use mod_precdefs
  use mod_params
  use time_vars, only: dt
  use field_vars, only: density_method => has_density, viscosity_method => has_viscosity
  use mod_physics, only: seawater_density, seawater_viscosity
  use mod_particle, only: t_particle
  use mod_fieldset, only: t_fieldset
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: vertical_velocity
  !---------------------------------------------
  ! [variable/type definition]
  !===================================================
contains
  !===========================================
  subroutine vertical_velocity(p, fieldset, time)
    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)            :: time
    real(rk)                        :: vert_vel

    dbghead(vertical_velocity)
    
    vert_vel = buoyancy(p, fieldset, time, p%delta_rho) + resuspension(p, fieldset, time)
    
    p%depth1 = p%depth1 + (vert_vel * dt)
    p%w1 = p%w1 + vert_vel
    p%vel_vertical = vert_vel
    
    debug(p%depth1)
    debug(p%w1)
    
    dbgtail(vertical_velocity)
    return
  end subroutine vertical_velocity
  !===========================================
  real(rk) function resuspension(p, fieldset, time) result(res)
    type(t_particle), intent(in) :: p
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in) :: time
    real(rk) :: i, j, k
    real(rk) :: u, v

    res = ZERO
    if (p%state /= BOTTOM) return

    if (resuspension_coeff >= ZERO) then

      i = p%ir0
      j = p%jr0
      k = p%kr0

      u = fieldset%get("U", time, i, j, k)
      v = fieldset%get("V", time, i, j, k)

      res = sqrt(u**2.+v**2.) * resuspension_coeff
    end if

    return
  end function resuspension
  !===========================================
  real(rk) function buoyancy(p, fieldset, time, delta_rho) result(res)
    !---------------------------------------------
    ! Calculate the vertical velocity due to buoyancy
    ! TODO: Which timestep should be used? (original or t + dt?)
    !---------------------------------------------
    type(t_particle), intent(in)    :: p
    type(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)            :: time
    real(rk)                        :: i, j, k
    real(rk)                        :: rho_sw ! Fluid density
    real(rk), intent(out)           :: delta_rho ! Density difference
    real(rk)                        :: kin_visc

    dbghead(buoyancy)

    debug(time)
    debug(p%depth0)
    debug(p%depth1)
    debug(p%w1)

    i = p%ir0
    j = p%jr0
    k = p%kr0

    rho_sw = ONE

    ! Density
    rho_sw = seawater_density(fieldset, time, i, j, k, density_method, p%depth0)
    delta_rho = p%rho - rho_sw

    ! Viscosity
    kin_visc = seawater_viscosity(fieldset, time, i, j, k, viscosity_method)

    res = Kooi_vertical_velocity(delta_rho, p%radius, rho_sw, kin_visc)

    dbgtail(buoyancy)
    return
  end function buoyancy
  !===========================================
  real(rk) function Kooi_vertical_velocity(delta_rho, rad_p, rho_env, kin_visc) result(res)
    !---------------------------------------------
    ! Calculate vertical velocity
    ! Reference: Kooi 2017
    !---------------------------------------------
    real(rk), intent(in) :: delta_rho, rad_p, rho_env
    real(rk), intent(in) :: kin_visc  ! Kinematic viscosity
    real(rk)             :: d_star    ! Dimensionless diameter
    real(rk)             :: w_star    ! Dimensionless settling velocity
    ! real(rk)             ::  ! Density difference

    dbghead(Kooi_vertical_velocity)

    debug(rad_p); 
    debug(rho_env)

    ! kin_visc = mu_env / rho_env ! NOT USING VISCOSITY???
    ! delta_rho = rho_p - rho_env

    debug(kin_visc); 
    debug(delta_rho)

    res = ZERO

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

end module mod_vertical_motion
