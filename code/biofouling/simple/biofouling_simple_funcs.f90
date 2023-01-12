#warning "This module has not been tested"
submodule(mod_biofouling) mod_biofouling_funcs
  !----------------------------------------------------------------
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
  module procedure get_ambient_chla_int_idx
  real(rk) :: lpp, spp, cya

  lpp = fieldset%get("LPP", time, i, j, k)
  spp = fieldset%get("SPP", time, i, j, k)
  cya = fieldset%get("CYA", time, i, j, k)

  res = 1.6 * 1.e6 * (spp + lpp + cya)

  end procedure get_ambient_chla_int_idx
  !===========================================
  module procedure get_ambient_chla_real_idx
  real(rk) :: lpp, spp, cya

  lpp = fieldset%get("LPP", time, i, j, k)
  spp = fieldset%get("SPP", time, i, j, k)
  cya = fieldset%get("CYA", time, i, j, k)

  res = 1.6 * 1.e6 * (spp + lpp + cya)

  end procedure get_ambient_chla_real_idx
end submodule mod_biofouling_funcs
