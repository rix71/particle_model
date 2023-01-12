#warning "This module has not been tested"
#include "cppdefs.h"
#include "file.h"
module mod_biofouling
  !----------------------------------------------------------------
  ! Biological growth on microplastic particles. Simple model that
  ! uses chl-a concentration as a proxy for the growth of biofilm.
  ! Reference:
  ! ! Add reference here
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  use run_params, only: biofouling_nmlfilename
  use mod_particle, only: t_particle
  use mod_fieldset, only: t_fieldset
  use time_vars, only: dt
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: init_biofouling, biofouling
  !---------------------------------------------
  real(rk) :: growth_init_threshold
  real(rk) :: h_bf_max
  real(rk) :: growth_timescale
  real(rk) :: rho_bf
  !---------------------------------------------
  interface
    module subroutine set_biofouling_fields(fieldset)
      type(t_fieldset), intent(inout) :: fieldset
    end subroutine set_biofouling_fields
    !---------------------------------------------
    module function get_ambient_chla_int_idx(fieldset, time, i, j, k) result(res)
      type(t_fieldset), intent(in) :: fieldset
      real(rk), intent(in) :: time
      integer, intent(in) :: i, j, k
      real(rk) :: res
    end function get_ambient_chla_int_idx
    !---------------------------------------------
    module function get_ambient_chla_real_idx(fieldset, time, i, j, k) result(res)
      type(t_fieldset), intent(in) :: fieldset
      real(rk), intent(in) :: time
      real(rk), intent(in) :: i, j, k
      real(rk) :: res
    end function get_ambient_chla_real_idx
  end interface
  !---------------------------------------------
  interface get_ambient_chla
    module procedure get_ambient_chla_int_idx
    module procedure get_ambient_chla_real_idx
  end interface get_ambient_chla
  !---------------------------------------------
  integer :: ierr
  !===================================================
contains
  !===========================================
  subroutine init_biofouling(fieldset)
    type(t_fieldset), intent(inout) :: fieldset
    namelist /biofouling/ growth_init_threshold, growth_timescale, h_bf_max, rho_bf

    FMT1, "======== Init biofouling ========"
    open (NMLFILE, file=trim(biofouling_nmlfilename), action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error('biofouling :: init_biofouling', "Failed to open "//trim(biofouling_nmlfilename), ierr)
    read (NMLFILE, nml=biofouling, iostat=ierr)
    if (ierr .ne. 0) call throw_error('biofouling :: init_biofouling', "Failed to read "//trim(biofouling_nmlfilename), ierr)
    close (NMLFILE)
    if (ierr .ne. 0) call throw_error('biofouling :: init_biofouling', "Failed to close "//trim(biofouling_nmlfilename), ierr)

    FMT2, LINE
    FMT2, "&biofouling"
    FMT3, var2val(growth_init_threshold)
    FMT3, var2val(growth_timescale)
    FMT3, var2val(h_bf_max)
    FMT3, var2val(rho_bf)

    call set_biofouling_fields(fieldset)

    return
  end subroutine init_biofouling
  !===========================================
  subroutine biofouling(p, fieldset, time)
    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in) :: time
    real(rk) :: chla, bf_growth

    chla = get_ambient_chla(fieldset, time, p%i0, p%j0, p%k0)
    ! If chl-a concentration is above a threshold, biofouling can occur
    if (chla > growth_init_threshold) then
      bf_growth = (h_bf_max - p%h_biofilm) / growth_timescale * dt
      p%h_biofilm = p%h_biofilm + bf_growth
      p%radius = p%radius + p%h_biofilm
      p%rho = (p%radius0**3.*p%rho0 + (p%radius**3.-p%radius0**3.) * rho_bf) / (p%radius)**3.
    end if

  end subroutine biofouling

end module mod_biofouling
