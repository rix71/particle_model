#include "cppdefs.h"
#include "particle.h"
module mod_vertical_motion
  !----------------------------------------------------------------
  ! Calculates the particles' vertical velocity
  !----------------------------------------------------------------
  use mod_errors
  use mod_precdefs
  use mod_params
  use time_vars, only: dt
  ! use field_vars, only: density_method, viscosity_method
  use mod_physics, only: seawater_density, seawater_viscosity, bottom_friction_velocity
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

    vert_vel = buoyancy(p, fieldset, time, p%delta_rho, p%kin_visc) + resuspension(p, fieldset, time, p%u_star)

    p%depth1 = p%depth1 + (vert_vel * dt)
    p%w1 = p%w1 + vert_vel
    p%vel_vertical = vert_vel

    return
  end subroutine vertical_velocity
  !===========================================
  real(rk) function resuspension(p, fieldset, time, u_star) result(res)
    !---------------------------------------------
    ! Gives a settled particle vertical velocity if
    ! the bottom friction velocity exceeds a certain threshold
    ! TODO: currently, resuspension_threshold is a namelist parameter,
    ! perhaps should be calculated as the particles' critical flow velocity.
    ! Ref: Erosion Behavior of Different Microplastic Particles in Comparison to Natural Sediments
    !       Kryss Waldschläger and Holger Schüttrumpf
    !       Environmental Science & Technology 2019 53 (22), 13219-13227
    !       DOI: 10.1021/acs.est.9b05394
    !---------------------------------------------
    type(t_particle), intent(in) :: p
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in) :: time
    real(rk) :: i, j, k
    real(rk), intent(out) :: u_star

    res = ZERO
    u_star = ZERO
    if (p%state /= ST_BOTTOM) then
      return
    end if

    if (resuspension_coeff >= ZERO) then

      i = p%ir0
      j = p%jr0
      k = p%kr0

      u_star = bottom_friction_velocity(fieldset, time, i, j, k)
      if (u_star > resuspension_threshold) then
        res = u_star * resuspension_coeff
      end if
    end if

    return
  end function resuspension
  !===========================================
  real(rk) function buoyancy(p, fieldset, time, delta_rho, kin_visc) result(res)
    !---------------------------------------------
    ! Calculate the vertical velocity due to buoyancy
    ! TODO: Which timestep should be used? (original or t + dt?)
    !---------------------------------------------
    type(t_particle), intent(in)    :: p
    type(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)            :: time
    real(rk)                        :: i, j, k
    real(rk)                        :: rho_sw    ! Fluid density
    real(rk), intent(out)           :: delta_rho ! Density difference
    real(rk), intent(out)           :: kin_visc  ! Kinematic viscosity (surrounding the particle)

    res = ZERO
    if (p%state /= ST_SUSPENDED) then
      return
    end if

    i = p%ir0
    j = p%jr0
    k = p%kr0

    rho_sw = ONE

    ! Density
    rho_sw = seawater_density(fieldset, time, i, j, k, p%depth0)
    delta_rho = p%rho - rho_sw

    ! Viscosity
    kin_visc = seawater_viscosity(fieldset, time, i, j, k, p%depth0)

    res = Kooi_vertical_velocity(delta_rho, p%radius, rho_sw, kin_visc)

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

    ! kin_visc = mu_env / rho_env ! NOT USING VISCOSITY???
    ! delta_rho = rho_p - rho_env

    dbghead(Kooi_vertical_velocity)
    debug(delta_rho)
    debug(kin_visc)

    res = ZERO
    d_star = (delta_rho * g * (2.*rad_p)**3.) / (rho_env * (kin_visc**2.)) ! g negative?
    debug(d_star)
    if (d_star < 0.05) then
      DBG, "d_star < 0.05"
      w_star = 1.74e-4 * (d_star**2)
    else if (d_star > 5.e9) then
      DBG, "d_star > 5.e9"
      w_star = 1000.
    else
      DBG, "d_star in range"
      w_star = 10.**(-3.76715 + (1.92944 * log10(d_star)) - (0.09815 * log10(d_star)**2.) &
                     - (0.00575 * log10(d_star)**3.) + (0.00056 * log10(d_star)**4.))
    end if
    debug(w_star)
    if (delta_rho > ZERO) then
      DBG, "delta_rho > 0"
      res = -1.0 * ((delta_rho / rho_env) * g * w_star * kin_visc)**(1./3.) ! Getting NaNs with -1*g
    else
      DBG, "delta_rho < 0"
      res = (-1.0 * (delta_rho / rho_env) * g * w_star * kin_visc)**(1./3.)
    end if
    debug(res)
    dbgtail(Kooi_vertical_velocity)
    return
  end function Kooi_vertical_velocity

end module mod_vertical_motion
