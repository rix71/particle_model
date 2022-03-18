#include "cppdefs.h"
module mod_advection
  use mod_precdefs
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
  subroutine advect_2d(p, fieldset, time)
    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)            :: time
    real(rk)                        :: x0, x1, x2, &
                                       y0, y1, y2, &
                                       lon1, lat1, &
                                       u1, u2, &
                                       v1, v2
    real(rk)                        :: i, j

    dbghead(advect_2d)

    ! x, y --> i, j
    i = p%ir0
    j = p%jr0
    ! call fieldset%search_indices(x=p%lon0, y=p%lat0, ir=i, jr=j)
    call fieldset%domain%lonlat2xy(p%lon0, p%lat0, x0, y0)

    debug(x0); debug(p%lon0); debug(i)
    debug(y0); debug(p%lat0); debug(j)

    u1 = fieldset%get("U", time, i, j)
    v1 = fieldset%get("V", time, i, j)

    x1 = x0 + u1 * dt
    y1 = y0 + v1 * dt

    debug(x1)
    debug(y1)

    ! x1, y1 --> i, j
    call fieldset%domain%xy2lonlat(x1, y1, lon1, lat1)
    call fieldset%search_indices(x=lon1, y=lat1, ir=i, jr=j)

    ! time --> time+1
    u2 = fieldset%get("U", time + dt, i, j)
    v2 = fieldset%get("V", time + dt, i, j)

    x2 = x0 + 0.5 * (u1 + u2) * dt
    y2 = y0 + 0.5 * (v1 + v2) * dt

    call fieldset%domain%xy2lonlat(x2, y2, p%lon1, p%lat1)
    call fieldset%search_indices(x=p%lon1, y=p%lat1, i=p%i1, j=p%j1, ir=p%ir1, jr=p%jr1)

    p%u1 = 0.5 * (u1 + u2)
    p%v1 = 0.5 * (v1 + v2)

    debug(x2); debug(p%lon1)
    debug(y2); debug(p%lat1)

    dbgtail(advect_2d)
    return
  end subroutine advect_2d
  !===========================================
  subroutine advect_3d(p, fieldset, time)
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

    dbghead(advect_3d)
    debug(p%lon0); debug(p%lat0); debug(p%depth0); 
    i = p%ir0
    j = p%jr0
    k = p%kr0
    ! call fieldset%search_indices(time, p%lon0, p%lat0, p%depth0, ir=i, jr=j, kr=k)
    debug(i); debug(j); debug(k); 
    call fieldset%domain%lonlat2xy(p%lon0, p%lat0, x0, y0)
    z0 = p%depth0

    debug(x0); debug(y0); debug(z0)

    u1 = fieldset%get("U", time, i=i, j=j, k=k)
    v1 = fieldset%get("V", time, i=i, j=j, k=k)
    w1 = fieldset%get("W", time, i=i, j=j, k=k)

    x1 = x0 + u1 * dt
    y1 = y0 + v1 * dt
    z1 = z0 + w1 * dt

    debug(x1); debug(y1); debug(z1)

    call fieldset%domain%xy2lonlat(x1, y1, lon1, lat1)
    call fieldset%search_indices(t=time + dt, x=lon1, y=lat1, z=z1, ir=i, jr=j, kr=k)

    u2 = fieldset%get("U", time + dt, i, j, k)
    v2 = fieldset%get("V", time + dt, i, j, k)
    w2 = fieldset%get("W", time + dt, i, j, k)

    x2 = x0 + 0.5 * (u1 + u2) * dt
    y2 = y0 + 0.5 * (v1 + v2) * dt
    z2 = z0 + 0.5 * (w1 + w2) * dt

    debug(x2); debug(y2); debug(z2)

    call fieldset%domain%xy2lonlat(x2, y2, p%lon1, p%lat1)
    call fieldset%search_indices(t=time + dt, x=p%lon1, y=p%lat1, z=z2, i=p%i1, j=p%j1, k=p%k1, ir=p%ir1, jr=p%jr1, kr=p%kr1)
    p%depth1 = z2

    p%u1 = 0.5 * (u1 + u2)
    p%v1 = 0.5 * (v1 + v2)
    p%w1 = 0.5 * (w1 + w2)

    dbgtail(advect_3d)
    return
  end subroutine advect_3d
  !===========================================
  subroutine advect(p, fieldset, time, adv_3d)
    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)            :: time
    logical, intent(in)             :: adv_3d

    select case (adv_3d)
    case (.true.)
      call advect_3d(p, fieldset, time)
    case (.false.)
      call advect_2d(p, fieldset, time)
    end select

    return
  end subroutine advect
  !===========================================

end module mod_advection
