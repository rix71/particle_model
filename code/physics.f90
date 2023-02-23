#ifdef SMAGORINSKY_FULL_FIELD
#error SMAGORINSKY_FULL_FIELD not implemented
#endif
#include "cppdefs.h"
#include "field.h"
module mod_physics
  !----------------------------------------------------------------
  ! This module contains functions to calculate physical parameters
  ! TODO:
  ! - drag
  !----------------------------------------------------------------
  use mod_errors
  use mod_interp
  use mod_precdefs
  use mod_params
  use field_vars, only: has_bottom_stress, density_method, viscosity_method, zax_style
  use time_vars, only: dt
  use mod_fieldset, only: t_fieldset
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: seawater_viscosity, seawater_density, diffusion_brown, &
            light_intensity, Ah_Smagorinsky, bottom_friction_velocity, normal_random
  !---------------------------------------------
  ! integer :: ierr ! Unused for now
  !===================================================
contains
  !===========================================
  real(rk) function seawater_viscosity(fieldset, time, i, j, k, depth) result(res)
    !---------------------------------------------
    ! Kinematic viscosity of seawater
    !---------------------------------------------
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in)         :: time
    real(rk), intent(in)         :: i, j, k
    real(rk), intent(in)         :: depth

    res = ZERO
    select case (viscosity_method)
    case (VISC_VARIABLE)
      res = fieldset%get("VISC", time, i, j, k)
      if (isnan(res) .or. (res <= ZERO)) res = kin_visc_default
    case (VISC_CALC)
      res = seawater_viscosity_from_temp_and_salt(fieldset%get("TEMP", time, i, j, k), &
                                                  fieldset%get("SALT", time, i, j, k))
      res = res / seawater_density(fieldset, time, i, j, k, depth)
    case (VISC_DEFAULT)
      res = kin_visc_default
    case default
      call throw_error("physics :: seawater_viscosity", "Viscosity method not recognized!")
    end select

    return
  end function seawater_viscosity
  !===========================================
  real(rk) function seawater_density(fieldset, time, i, j, k, depth) result(res)
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in)         :: time
    real(rk), intent(in)         :: i, j, k
    real(rk), intent(in)         :: depth

    res = ZERO
    select case (density_method)
    case (RHO_VARIABLE)
      res = fieldset%get("RHO", time, i, j, k)
    case (RHO_CALC)
      res = seawater_density_from_temp_and_salt(fieldset%get("TEMP", time, i, j, k), &
                                                fieldset%get("SALT", time, i, j, k), &
                                                depth)
    case (RHO_DEFAULT)
      res = sw_rho_default
    case default
      call throw_error("physics :: seawater_density", "Density method not recognized!")
    end select

    return
  end function seawater_density
  !===========================================
  real(rk) function diffusion_brown(T, r, mu) result(res)
    real(rk), intent(in) :: T, r, mu

    res = (k_b * (T + 273.16)) / (6.*pi * mu * r)

    return
  end function diffusion_brown
  !===========================================
  real(rk) function light_intensity(t, z) result(res)
    real(rk), intent(in) :: t ! Hour of day
    real(rk), intent(in) :: z ! Depth
    real(rk) :: I0 ! Light intensity at sea surface

    I0 = Im * sin(2.*pi * (t / 24.-(6./24.)))
    if (I0 < ZERO) I0 = ZERO
    res = I0 * exp(eps_light * z)

    return
  end function light_intensity
  !===========================================
  function Ah_Smagorinsky(fieldset, time, i, j, k) result(Ah_s)
    type(t_fieldset), intent(in)  :: fieldset
    real(rk), intent(in)          :: time
    real(rk), intent(in)           :: i, j
    real(rk), intent(in), optional :: k
    real(rk)                      :: Ah_s
    real(rk)                      :: dx, dy
    real(rk)                      :: i0, j0, i1, j1, ih, jh, kh
    real(rk)                      :: uh0, vh0, uh1, vh1, u0h, v0h, u1h, v1h
    real(rk)                      :: dudx, dudy, dvdx, dvdy

    i0 = real(i, rk); 
    j0 = real(j, rk); 
    i1 = real(i, rk) + ONE; 
    j1 = real(j, rk) + ONE; 
    ih = real(i, rk) + HALF; 
    jh = real(j, rk) + HALF; 
    if (present(k)) then

      kh = real(k, rk) + HALF; 
      uh0 = fieldset%get("U", time, ih, j0, kh); 
      vh0 = fieldset%get("V", time, ih, j0, kh); 
      uh1 = fieldset%get("U", time, ih, j1, kh); 
      vh1 = fieldset%get("V", time, ih, j1, kh); 
      u0h = fieldset%get("U", time, i0, jh, kh); 
      v0h = fieldset%get("V", time, i0, jh, kh); 
      u1h = fieldset%get("U", time, i1, jh, kh); 
      v1h = fieldset%get("V", time, i1, jh, kh); 
    else
      uh0 = fieldset%get("U", time, ih, j0); 
      vh0 = fieldset%get("V", time, ih, j0); 
      uh1 = fieldset%get("U", time, ih, j1); 
      vh1 = fieldset%get("V", time, ih, j1); 
      u0h = fieldset%get("U", time, i0, jh); 
      v0h = fieldset%get("V", time, i0, jh); 
      u1h = fieldset%get("U", time, i1, jh); 
      v1h = fieldset%get("V", time, i1, jh); 
    end if

    dx = fieldset%domain%dx
    dy = fieldset%domain%dy

    dudx = (u1h - u0h) / dx; 
    dudy = (uh1 - uh0) / dy; 
    dvdx = (v1h - v0h) / dx; 
    dvdy = (vh1 - vh0) / dy; 
    ! Ah_s = Cm_smagorinsky * sqrt(HALF * (dudx - dvdy)**2 + HALF * (dudy + dvdx)**2)
    Ah_s = Cm_smagorinsky * (dx * dy) * sqrt(dudx**2 + dvdy**2 + HALF * (dudy + dvdx)**2.)

    return
  end function Ah_Smagorinsky
  !===========================================
  subroutine bottom_stress(fieldset, time, i, j, k, taubx, tauby)
    !---------------------------------------------
    ! Calculate the bottom stress using a logarithmic profile
    ! Returns taubx/rho, tauby/rho (stress divided by density)
    !! Still debugging this
    !---------------------------------------------
    type(t_fieldset), intent(in)  :: fieldset
    real(rk), intent(in)          :: time
    real(rk), intent(in)          :: i, j
    real(rk), intent(in)          :: k
    real(rk), intent(out)         :: taubx, tauby
    real(rk)                      :: u, v, h
    real(rk)                      :: C ! Drag coefficient
    real(rk), parameter           :: von_Karman = 0.4_rk
    real(rk), parameter           :: c_drag_min = 0.0025_rk ! Minimum drag coefficient
    real(rk), allocatable         :: zax(:)
    real(rk)                      :: level_idx

    ! if (int(k) /= fieldset%zax_bot_idx) call throw_error("physics :: bottom_stress", "Not bottom layer index!")

    level_idx = real(int(k), rk)
    if (fieldset%bottom_is_nan("U")) level_idx = level_idx + ONE

    if (level_idx < fieldset%zax_bot_idx .or. level_idx > fieldset%zax_top_idx) then
      ERROR, "physics :: bottom_stress :: level_idx out of range"
      ERROR, "level_idx = ", level_idx
      ERROR, "bottom_nan = ", fieldset%bottom_is_nan("U")
    end if

    u = fieldset%get("U", time, i, j, k=level_idx)
    v = fieldset%get("V", time, i, j, k=level_idx)
    ! Automatic allocation should happen here
    select case (zax_style)
    case (STATIC_DEPTH_VALUES)
      zax = fieldset%get_zax()
    case (DEPTH_VALUES, LAYER_THICKNESS)
      zax = fieldset%get_zax(time, int(i), int(j))
    end select
    ! h is the height of the layer above the bottom
    h = -zax(int(level_idx))
    C = max((von_Karman**2./((log(h / roughness_height))**2.)), c_drag_min)
    taubx = -C * sqrt(u**2.+v**2.) * u
    tauby = -C * sqrt(u**2.+v**2.) * v

    return
  end subroutine bottom_stress
  !===========================================
  real(rk) function bottom_friction_velocity(fieldset, time, i, j, k) result(u_star)
    !---------------------------------------------
    ! Calculate the bottom friction velocity u_star from the velocity components
    ! of the above layer
    !---------------------------------------------
    type(t_fieldset), intent(in)  :: fieldset
    real(rk), intent(in)          :: time
    real(rk), intent(in)          :: i, j
    real(rk), intent(in)          :: k
    real(rk)                      :: taubx, tauby ! Bottom stress

    if (has_bottom_stress) then
      taubx = fieldset%get("TAUBX", time, i, j)
      tauby = fieldset%get("TAUBY", time, i, j)
    else
      call bottom_stress(fieldset, time, i, j, k, taubx, tauby)
    end if
    u_star = (taubx**2.+tauby**2.)**0.25

    return
  end function bottom_friction_velocity
  !===========================================
  real(rk) function seawater_viscosity_from_temp_and_salt(T, S)
    !---------------------------------------------
    ! Dynamic viscosity of seawater
    ! TODO: Add proper reference
    ! https://ittc.info/media/4048/75-02-01-03.pdf
    !---------------------------------------------
    real(rk), intent(in) :: T, S
    real(rk) :: A, B, mu

    A = 1.541_rk + 1.998e-2_rk * T - 9.52e-5_rk * T**2.
    B = 7.974_rk - 7.561e-2_rk * T + 4.724e-4_rk * T**2.
    mu = 4.2844e-5_rk + (ONE / (0.156_rk * (T + 64.993_rk)**2.-91.296_rk))
    seawater_viscosity_from_temp_and_salt = mu * (ONE + A * S + B * S**2.)

    return
  end function seawater_viscosity_from_temp_and_salt
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
