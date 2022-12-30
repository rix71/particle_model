#include "cppdefs.h"
#include "particle.h"
#include "advection.h"
module mod_advection
  use mod_precdefs
  use mod_errors
  use time_vars, only: dt
  use mod_particle, only: t_particle
  use mod_fieldset, only: t_fieldset
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: advect
  !---------------------------------------------
  !===================================================
contains
  !===========================================
  subroutine advect_EE_2D(p, fieldset, time)
    !---------------------------------------------
    ! Explicit Euler scheme for advection (2D)
    !---------------------------------------------
    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)            :: time
    real(rk)                        :: x0, x1, &
                                       y0, y1, &
                                       u, v
    real(rk)                        :: i, j

    i = p%ir0
    j = p%jr0
    call fieldset%domain%lonlat2xy(p%lon0, p%lat0, x0, y0)

    u = fieldset%get("U", time, i, j)
    v = fieldset%get("V", time, i, j)

    x1 = x0 + u * dt
    y1 = y0 + v * dt

    call fieldset%domain%xy2lonlat(x1, y1, p%lon1, p%lat1)
    call fieldset%search_indices(x=p%lon1, y=p%lat1, i=p%i1, j=p%j1, ir=p%ir1, jr=p%jr1)

    p%u1 = u
    p%v1 = v

    return
  end subroutine advect_EE_2D
  !===========================================
  subroutine advect_EE_3D(p, fieldset, time)
    !---------------------------------------------
    ! Explicit Euler scheme for advection (3D)
    !---------------------------------------------
    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)            :: time
    real(rk)                        :: x0, x1, &
                                       y0, y1, &
                                       z0, z1, &
                                       u, v, w
    real(rk)                        :: i, j, k

    i = p%ir0
    j = p%jr0
    k = p%kr0
    call fieldset%domain%lonlat2xy(p%lon0, p%lat0, x0, y0)
    z0 = p%depth0

    u = fieldset%get("U", time, i, j, k)
    v = fieldset%get("V", time, i, j, k)
#ifndef NO_ADVECT_VERTICAL
    w = fieldset%get("W", time, i, j, k)
#else
    w = ZERO
#endif

    x1 = x0 + u * dt
    y1 = y0 + v * dt
    z1 = z0 + w * dt

    call fieldset%domain%xy2lonlat(x1, y1, p%lon1, p%lat1)
    call fieldset%search_indices(t=time + dt, x=p%lon1, y=p%lat1, i=p%i1, j=p%j1, k=p%k1, ir=p%ir1, jr=p%jr1, kr=p%kr1)

    p%u1 = u
    p%v1 = v
    p%w1 = w

    return
  end subroutine advect_EE_3D
  !===========================================
  subroutine advect_RK2_2D(p, fieldset, time)
    !---------------------------------------------
    ! Second order Runge-Kutta scheme for advection (2D)
    !---------------------------------------------
    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)            :: time
    real(rk)                        :: x0, x1, x2, &
                                       y0, y1, y2, &
                                       lon1, lat1, &
                                       u1, u2, &
                                       v1, v2
    real(rk)                        :: i, j

    i = p%ir0
    j = p%jr0
    call fieldset%domain%lonlat2xy(p%lon0, p%lat0, x0, y0)

    u1 = fieldset%get("U", time, i, j)
    v1 = fieldset%get("V", time, i, j)

    x1 = x0 + u1 * dt
    y1 = y0 + v1 * dt

    call fieldset%domain%xy2lonlat(x1, y1, lon1, lat1)
    call fieldset%search_indices(x=lon1, y=lat1, ir=i, jr=j)

    u2 = fieldset%get("U", time + dt, i, j)
    v2 = fieldset%get("V", time + dt, i, j)

    x2 = x0 + 0.5 * (u1 + u2) * dt
    y2 = y0 + 0.5 * (v1 + v2) * dt

    call fieldset%domain%xy2lonlat(x2, y2, p%lon1, p%lat1)
    call fieldset%search_indices(x=p%lon1, y=p%lat1, i=p%i1, j=p%j1, ir=p%ir1, jr=p%jr1)

    p%u1 = 0.5 * (u1 + u2)
    p%v1 = 0.5 * (v1 + v2)

    return
  end subroutine advect_RK2_2D
  !===========================================
  subroutine advect_RK2_3D(p, fieldset, time)
    !---------------------------------------------
    ! Second order Runge-Kutta scheme for advection (3D)
    !---------------------------------------------
    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)            :: time
    real(rk)                        :: x0, x1, x2, &
                                       y0, y1, y2, &
                                       z0, z1, z2, &
                                       lon1, lat1, &
                                       u1, u2, &
                                       v1, v2, &
                                       w1, w2
    real(rk)                        :: i, j, k

    i = p%ir0
    j = p%jr0
    k = p%kr0

    call fieldset%domain%lonlat2xy(p%lon0, p%lat0, x0, y0)
    z0 = p%depth0

    u1 = fieldset%get("U", time, i=i, j=j, k=k); if (isnan(u1)) u1 = ZERO
    v1 = fieldset%get("V", time, i=i, j=j, k=k); if (isnan(v1)) v1 = ZERO
#ifndef NO_ADVECT_VERTICAL
    w1 = fieldset%get("W", time, i=i, j=j, k=k); if (isnan(w1)) w1 = ZERO
#else
    w1 = ZERO
#endif

    x1 = x0 + u1 * dt
    y1 = y0 + v1 * dt
    z1 = z0 + w1 * dt

    call fieldset%domain%xy2lonlat(x1, y1, lon1, lat1)
    call fieldset%search_indices(t=time + dt, x=lon1, y=lat1, z=z1, ir=i, jr=j, kr=k)

    u2 = fieldset%get("U", time + dt, i, j, k); if (isnan(u2)) u2 = ZERO
    v2 = fieldset%get("V", time + dt, i, j, k); if (isnan(v2)) v2 = ZERO
#ifndef NO_ADVECT_VERTICAL
    w2 = fieldset%get("W", time + dt, i, j, k); if (isnan(w2)) w2 = ZERO
#else
    w2 = ZERO
#endif

    x2 = x0 + 0.5 * (u1 + u2) * dt
    y2 = y0 + 0.5 * (v1 + v2) * dt
    z2 = z0 + 0.5 * (w1 + w2) * dt

    call fieldset%domain%xy2lonlat(x2, y2, p%lon1, p%lat1)
    call fieldset%search_indices(t=time + dt, x=p%lon1, y=p%lat1, z=z2, i=p%i1, j=p%j1, k=p%k1, ir=p%ir1, jr=p%jr1, kr=p%kr1)
    p%depth1 = z2

    p%u1 = 0.5 * (u1 + u2)
    p%v1 = 0.5 * (v1 + v2)
    p%w1 = 0.5 * (w1 + w2)

    return
  end subroutine advect_RK2_3D
  !===========================================
  subroutine advect(p, fieldset, time, adv_method, adv_3d)
    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)            :: time
    integer, intent(in)             :: adv_method
    logical, intent(in)             :: adv_3d

    if (p%state /= ST_SUSPENDED) return

    select case (adv_method)
    case (ADV_NONE)
      return
    case (ADV_RK2)
      select case (adv_3d)
      case (.true.)
        call advect_RK2_3D(p, fieldset, time)
      case (.false.)
        call advect_RK2_2D(p, fieldset, time)
      end select
    case (ADV_EE)
      select case (adv_3d)
      case (.true.)
        call advect_EE_3D(p, fieldset, time)
      case (.false.)
        call advect_EE_2D(p, fieldset, time)
      end select
    case default
      call throw_error("advection :: advect", "Undefined advection method!")
    end select

    return
  end subroutine advect

end module mod_advection
