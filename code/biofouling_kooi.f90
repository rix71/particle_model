#warning This module is not tested
#include "cppdefs.h"
#include "file.h"
module mod_biofouling
  !----------------------------------------------------------------
  ! Biological growth on microplastic particles. Causes density changes.
  ! References:
  ! [1] Global modeled sinking characteristics of biofouled microplastic. Journal of Geophysical Research:
  !     Lobelle, D., Kooi, M., Koelmans, A. A., Laufkötter, C., Jongedijk, C. E., Kehl, C., & van Sebille, E. (2021).
  !     Oceans, 126, e2020JC017098. https://doi.org/10.1029/2020JC017098
  !
  ! [2] Ups and Downs in the Ocean: Effects of Biofouling on Vertical Transport of Microplastics
  !     Merel Kooi, Egbert H. van Nes, Marten Scheffer, and Albert A. Koelmans
  !     Environmental Science & Technology 2017 51 (14), 7963-7971
  !     DOI: 10.1021/acs.est.6b04702
  !
  ! [3] Validation of a simple model accounting for light and temperature effect on microalgal growth.
  !     Bernard O, Rémond B.
  !     Bioresour Technol. 2012 Nov;123:520-7. doi: 10.1016/j.biortech.2012.07.022. Epub 2012 Jul 16. PMID: 22940363.
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  use mod_params, only: pi
  use mod_particle, only: t_particle
  use mod_fieldset, only: t_fieldset
  use mod_datetime, only: t_datetime
  use mod_physics, only: seawater_viscosity, seawater_density, diffusion_brown, light_intensity
  use time_vars, only: dt
  use field_vars, only: viscosity_method => has_viscosity, density_method => has_density
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: init_biofouling, biofouling
  !---------------------------------------------
  real(rk) :: r20, &              ! Respiration rate
              q10, &              ! Respiration increase per 10 deg
              vol_algal, &        ! Volume of algal cell
              gamma, &            ! Shear
              rho_bf, &           ! Density of biofilm
              growth_rate_max, &  ! Maximum growth rate under optimal conditions
              mortality_rate, &   ! Grazing mortality rate
              slope_init, &            ! Initial slope of growth rate
              I_opt, &            ! Optimal light intensity for algae growth
              T_min, T_max, T_opt ! Minimum, maximum and optimal temperature to sustain algae growth
  !---------------------------------------------
  interface
    module subroutine set_biofouling_fields(fieldset, filename)
      type(t_fieldset), intent(inout)       :: fieldset
      character(len=LEN_CHAR_L), intent(in) :: filename
    end subroutine set_biofouling_fields
    !---------------------------------------------
    module function get_alg_ambient(fieldset, time, i, j, k) result(res)
      type(t_fieldset), intent(in) :: fieldset
      real(rk), intent(in) :: time
      real(rk), intent(in) :: i, j, k
      real(rk)                     :: res
    end function get_alg_ambient
  end interface
  !---------------------------------------------
  integer :: ierr
  !===================================================
contains
  !===========================================
  subroutine init_biofouling(fieldset, filename)
    type(t_fieldset), intent(inout)       :: fieldset
    character(len=LEN_CHAR_L), intent(in) :: filename
    namelist /biofouling/ vol_algal, rho_bf, gamma, &
      growth_rate_max, q10, r20, mortality_rate, &
      slope_init, I_opt, T_min, T_max, T_opt

    FMT1, "======== Init biofouling ========"

    open (NMLFILE, file=NMLFILENAME, action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error('biofouling :: init_biofouling', "Failed to open "//NMLFILENAME, ierr)
    read (NMLFILE, nml=biofouling)
    close (NMLFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error('biofouling :: init_biofouling', "Failed to close "//NMLFILENAME, ierr)

    FMT2, LINE
    FMT2, "&biofouling"
    FMT3, var2val(vol_algal)
    FMT3, var2val(rho_bf)
    FMT3, var2val(gamma)
    FMT3, var2val(growth_rate_max)
    FMT3, var2val(q10)
    FMT3, var2val(r20)
    FMT3, var2val(mortality_rate)
    FMT3, var2val(slope_init)
    FMT3, var2val(I_opt)
    FMT3, var2val(T_min)
    FMT3, var2val(T_max)
    FMT3, var2val(T_opt)

    call set_biofouling_fields(fieldset, filename)

    return
  end subroutine init_biofouling
  !===========================================
  subroutine biofouling(p, fieldset, time, date)
    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)            :: time
    type(t_datetime), intent(in)    :: date
    real(rk)                        :: mu ! Dynamic viscosity
    real(rk)                        :: T  ! Temperature
    real(rk)                        :: a_coll, a_growth, a_mort, a_resp
    real(rk)                        :: v_bf, v_tot ! Volume of biofilm and total volume
    real(rk)                        :: alg_ambient ! Ambient algae concentration [1/m3]

    T = fieldset%get("TEMP", time, p%ir0, p%jr0, p%kr0)

    ! TODO: Consider making a function for this (seawater_dyn_viscosity)
    mu = seawater_viscosity(fieldset, time, p%ir0, p%jr0, p%kr0, viscosity_method) / &
         seawater_density(fieldset, time, p%ir0, p%jr0, p%kr0, density_method, p%depth0)

    alg_ambient = get_alg_ambient(fieldset, time, p%ir0, p%jr0, p%kr0)

    a_coll = (encounter_rate(p, T, mu) * alg_ambient) / p%surface_area()
    a_growth = primary_prod(T, p%depth0, real(date%hour, rk)) * p%aa_growth
    a_mort = mortality_rate * p%aa_growth
    a_resp = (q10**((T - 20.) / 10.)) * r20 * p%aa_growth

    p%aa_growth = p%aa_growth + (a_coll + a_growth - a_mort - a_resp) * dt

    v_bf = (vol_algal * p%aa_growth) * p%surface_area()
    v_tot = v_bf + p%volume()
    p%h_biofilm = ((v_tot * (3./(4.*pi)))**(1./3.)) - p%radius0

    p%radius = p%radius0 + p%h_biofilm

    p%rho = (p%radius0**3.*p%rho0 + (p%radius**3.-p%radius0**3.) * rho_bf) / (p%radius)**3.

    return
  end subroutine biofouling
  !===========================================
  real(rk) function primary_prod(T, depth, hour) result(res)
    !---------------------------------------------
    ! Primary productivity. Currently parametrised based on [3].
    ! Function of temperature and light intensity:
    !   mu(Tz, Iz) = mu_opt(Iz) * Phi(Tz),
    ! where mu_opt is the optimal growth rate under light conditions
    ! at depth z, and Phi is the temperature influence on the growth rate.
    !
    ! Ideally, this should only be an alternative if such data is not available from models.
    !---------------------------------------------
    real(rk), intent(in) :: T     ! Temperature at particle depth
    real(rk), intent(in) :: depth ! Depth of the particle
    real(rk), intent(in) :: hour  ! Hour of day
    real(rk)             :: Iz, mu_opt, phi

    Iz = light_intensity(hour, depth)
    mu_opt = (growth_rate_max * Iz) / (Iz + (growth_rate_max / slope_init) * (Iz / I_opt - 1)**2.)
    phi = ((T - T_max)*(T - T_min)**2.)/((T_opt - T_min) * ((T_opt - T_min)*(T - T_opt) - (T_opt - T_max)*(T_opt + T_min - 2.*T)))

    res = mu_opt * phi

  end function primary_prod
  !===========================================
  real(rk) function encounter_rate(p, T, mu) result(res)
    !---------------------------------------------
    ! The encounter kernel rate (beta_A) is the sum of Brownian motion (beta_Abrown),
    ! differential settling (beta_Asettling) and advective shear (beta_Ashear) collision frequencies [m3 s-1].
    !---------------------------------------------
    type(t_particle), intent(in) :: p
    real(rk), intent(in)         :: T  ! Temperature at particle depth
    real(rk), intent(in)         :: mu ! Dynamic viscosity
    real(rk)                     :: beta_Abrown, beta_Asettling, beta_Ashear
    real(rk)                     :: r_algae ! Radius of algal cell

    r_algae = ((3./4.) * (vol_algal / pi))**(1./3.)

    beta_Abrown = 4.*pi * (diffusion_brown(T, p%radius, mu) + diffusion_brown(T, r_algae, mu)) * (p%radius + r_algae)
    beta_Asettling = HALF * pi * p%radius**2 * abs(p%vel_vertical)
    beta_Ashear = 1.3 * gamma * (p%radius + r_algae)**3.

    res = beta_Abrown + beta_Asettling + beta_Ashear

    return
  end function encounter_rate
end module mod_biofouling
