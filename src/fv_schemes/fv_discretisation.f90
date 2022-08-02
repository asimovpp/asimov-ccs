!>  Submodule file fv_discretisation.smod
!
!>  @build discretisation
!
!>  Implementations of the finite volume method using the various discretisation schemes scheme

submodule(fv) fv_discretisation

  implicit none

contains

  !>  Calculates advection coefficient for neighbouring cell using CDS discretisation
  !
  !> @param[in] phi         - scalar field
  !> @param[in] mf          - mass flux at the face
  !> @param[in] bc          - flag indicating whether cell is on boundary
  !> @param[out] coeff      - advection coefficient to be calculated
  module subroutine calc_advection_coeff_cds(phi, mf, bc, coeff)
    type(central_field), intent(in) :: phi
    real(ccs_real), intent(in) :: mf
    integer(ccs_int), intent(in) :: bc
    real(ccs_real), intent(out) :: coeff

    real(ccs_real) :: interpolation_factor

    ! Dummy usage to prevent unused argument.
    associate (scalar => phi)
    end associate
    associate (mflux => mf)
    end associate

    if (bc == 0) then
      interpolation_factor = 0.5_ccs_real
    else
      interpolation_factor = 1.0_ccs_real
    end if
    coeff = interpolation_factor
  end subroutine calc_advection_coeff_cds

  !>  Calculates advection coefficient for neighbouring cell using UDS discretisation
  !
  !> @param[in] phi         - scalar field
  !> @param[in] mf          - mass flux at the face
  !> @param[in] bc          - flag indicating whether cell is on boundary
  !> @param[out] coeff      - advection coefficient to be calculated
  module subroutine calc_advection_coeff_uds(phi, mf, bc, coeff)
    type(upwind_field), intent(in) :: phi
    real(ccs_real), intent(in) :: mf
    integer(ccs_int), intent(in) :: bc
    real(ccs_real), intent(out) :: coeff

    ! Dummy usage to prevent unused argument.
    associate (scalar => phi, foo => bc)
    end associate

    if (mf < 0.0) then
      coeff = 1.0_ccs_real
    else
      coeff = 0.0_ccs_real
    end if

  end subroutine calc_advection_coeff_uds

end submodule fv_discretisation
