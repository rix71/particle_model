submodule(mod_biofouling) mod_biofouling_funcs
  !----------------------------------------------------------------
  ! The subroutine set_biofouling_fields should add necessary field/fields to
  ! the fieldset.
  !
  ! The function get_alg_ambient should handle retrieveing the value(s) at the particle location
  ! and calculations/conversions if necessary (return value must have the units [1/m-3]).
  !----------------------------------------------------------------
  implicit none
  !===================================================
  !===================================================
contains
  !===========================================
  module procedure set_biofouling_fields

  call fieldset%add_field("LPP", "iow_ergom_t_lpp")
  call fieldset%add_field("SPP", "iow_ergom_t_spp")
  call fieldset%add_field("CYA", "iow_ergom_t_cya")

  end procedure set_biofouling_fields
  !===========================================
  module procedure get_alg_ambient
  real(rk) :: lpp, spp, cya

  lpp = fieldset%get("LPP", time, i, j, k)
  spp = fieldset%get("SPP", time, i, j, k)
  cya = fieldset%get("CYA", time, i, j, k)

  res = 1.6 * 1.e6 * (spp + lpp + cya)

  end procedure get_alg_ambient
end submodule mod_biofouling_funcs
