#include "cppdefs.h"
module mod_biofouling
  !----------------------------------------------------------------
  ! Biological growth on microplastic particles. Causes density changes.
  ! Reference:
  !  Modeling the Pathways and Accumulation Patterns of Micro- and Macro-Plastics in the Mediterranean.
  !   Tsiaras K, Hatzonikolakis Y, Kalaroni S, Pollani A and Triantafyllou G (2021)
  !   Front. Mar. Sci. 8:743117. doi: 10.3389/fmars.2021.743117
  !----------------------------------------------------------------
  use mod_errors
  use mod_precdefs
  use mod_params, only: pi
  use time_vars, only: dt
  use mod_particle, only: t_particle
  use mod_fieldset, only: t_fieldset
  use nc_manager, only: nc_var_exists
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: init_biofouling, biofouling
  !---------------------------------------------
  real(rk) :: cell_size
  real(rk) :: Ba_mean
  real(rk) :: rho_biofilm
  real(rk) :: h_crit
  real(rk) :: D
  !---------------------------------------------
  character(LEN_CHAR_S) :: bavarname
  !---------------------------------------------
  integer :: ierr
  !===================================================
contains
  !===========================================
  subroutine init_biofouling(fieldset, filename)
    type(t_fieldset), intent(inout) :: fieldset
    character(len=LEN_CHAR_L) :: filename

    namelist /biofouling/ cell_size, Ba_mean, rho_biofilm, h_crit, D, bavarname

    FMT1, "======== Init biofouling ========"

    open (NMLFILE, file=NMLFILENAME, action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error('biofouling :: init_biofouling', "Failed to open "//NMLFILENAME, ierr)
    read (NMLFILE, nml=biofouling)
    close (NMLFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error('biofouling :: init_biofouling', "Failed to close "//NMLFILENAME, ierr)

    FMT2, LINE
    FMT2, "&biofouling"
    FMT3, var2val(cell_size)
    FMT3, var2val(Ba_mean)
    FMT3, var2val(rho_biofilm)
    FMT3, var2val(h_crit)
    FMT3, var2val(D)
    FMT3, var2val_char(bavarname)

    if (nc_var_exists(trim(filename), trim(bavarname))) then
      call fieldset%add_field("BA", trim(bavarname))
    else
      call throw_error("biofouling :: init_biofouling", "Variable for BA does not exist: "//trim(bavarname))
    end if

    return
  end subroutine init_biofouling
  !===========================================
  real(rk) function dh0(r) result(res)
    !---------------------------------------------
    ! Mean biofilm growth rate
    ! Parameters (from namelist, probably):
    !    D - diffusion coefficient for bacteria
    !    Ba_mean - mean cell abundance
    !    cell_size - cell size, typical value ~1 micrometer
    !---------------------------------------------
    real(rk), intent(in) :: r

    res = D / r * Ba_mean * pi / 4.*cell_size**3

    return
  end function dh0
  !===========================================
  real(rk) function detachment(p) result(res)
    !---------------------------------------------
    ! Detachment rate
    ! h_crit should be calculated
    !---------------------------------------------
    type(t_particle), intent(in) :: p

    res = max((0.92 * dh0(p%radius)) / (h_crit + 3.e-6), 0.025)

    return
  end function detachment
  !===========================================
  subroutine biofouling(p, fieldset, time)
    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in) :: time
    real(rk) :: i, j, k
    real(rk) :: Ba, &
                ! h, &
                growth

    i = p%ir0
    j = p%jr0
    k = p%kr0

    Ba = fieldset%get("BA", time, i, j, k)

    growth = ((dh0(p%radius) * Ba) / Ba_mean - detachment(p) * p%h_biofilm) * dt
    p%radius = max(p%radius + growth, p%radius0)
    p%h_biofilm = max(p%h_biofilm + growth, ZERO)
    p%rho = max(p%rho * (p%radius0**3 / p%radius**3) + rho_biofilm * (1 - (p%radius0**3 / p%radius**3)), p%rho0)

    return
  end subroutine biofouling

end module mod_biofouling
