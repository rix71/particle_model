#ifdef SMAGORINSKY_FULL_FIELD
#error SMAGORINSKY_FULL_FIELD not implemented
#endif
#include "cppdefs.h"
#include "particle.h"
module mod_diffusion
  use mod_precdefs
  use mod_errors
  use mod_params
  use time_vars, only: dt
  use mod_particle, only: t_particle
  use mod_fieldset, only: t_fieldset
  use mod_physics, only: normal_random, Ah_Smagorinsky
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: diffuse
  !---------------------------------------------
  !===================================================
contains
!===========================================
  subroutine diffuse_2D(p, fieldset, time)

    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)            :: time
    real(rk)                        :: Ah
    real(rk)                        :: i, j
    real(rk)                        :: x0, x1, &
                                       y0, y1

    i = p%ir1
    j = p%jr1

    call fieldset%domain%lonlat2xy(p%lon1, p%lat1, x0, y0)

#if defined(SMAGORINSKY_INTERP_UV)
    Ah = Ah_Smagorinsky(fieldset, time, i, j)
#elif defined(SMAGORINSKY_FULL_FIELD)
    Ah = get_Ah_Smagorinsky_full_field(fieldset, i, j)
#endif

    x1 = x0 + normal_random() * sqrt(2 * Ah * dt)
    y1 = y0 + normal_random() * sqrt(2 * Ah * dt)

    call fieldset%domain%xy2lonlat(x1, y1, p%lon1, p%lat1)
    call fieldset%search_indices(x=p%lon1, y=p%lat1, i=p%i1, j=p%j1, ir=p%ir1, jr=p%jr1)

    return
  end subroutine diffuse_2D
  !===========================================
  subroutine diffuse_3D(p, fieldset, time)

    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)            :: time
    real(rk)                        :: Ah, kv
    real(rk)                        :: i, j, k
    real(rk)                        :: x0, x1, &
                                       y0, y1, &
                                       z0, z1

    i = p%ir1
    j = p%jr1
    k = p%kr1

    call fieldset%domain%lonlat2xy(p%lon1, p%lat1, x0, y0)
    z0 = p%depth1

#if defined(SMAGORINSKY_INTERP_UV)
    Ah = Ah_Smagorinsky(fieldset, time, i, j, k)
#elif defined(SMAGORINSKY_FULL_FIELD)
    Ah = get_Ah_Smagorinsky_full_field(fieldset, i, j, k)
#endif
    kv = diffusion_vert_const

    x1 = x0 + normal_random() * sqrt(2 * Ah * dt)
    y1 = y0 + normal_random() * sqrt(2 * Ah * dt)
#ifndef NO_DIFFUSE_VERTICAL
    z1 = z0 + normal_random() * sqrt(2 * kv * dt)
#else
    z1 = z0
#endif

    call fieldset%domain%xy2lonlat(x1, y1, p%lon1, p%lat1)
    call fieldset%search_indices(x=p%lon1, y=p%lat1, i=p%i1, j=p%j1, ir=p%ir1, jr=p%jr1)

    return
  end subroutine diffuse_3D
  !===========================================
  subroutine diffuse(p, fieldset, time, dif_3d)

    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)            :: time
    logical, intent(in)             :: dif_3d

    if (p%state /= ST_SUSPENDED) return

    select case (dif_3d)
    case (.true.)
      call diffuse_3D(p, fieldset, time)
    case (.false.)
      call diffuse_2D(p, fieldset, time)
    end select

    return
  end subroutine diffuse
end module mod_diffusion
