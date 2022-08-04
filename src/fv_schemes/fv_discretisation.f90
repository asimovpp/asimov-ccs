!v Submodule file fv_discretisation.smod
!
!  Implementations of the finite volume method using the various discretisation schemes scheme
!
!  @build discretisation

submodule(fv) fv_discretisation

  implicit none

contains

  !> Calculates advection coefficient for neighbouring cell using CDS discretisation
  module subroutine calc_advection_coeff_cds(phi, mf, bc, coeff)
    type(central_field), intent(in) :: phi !< scalar field
    real(ccs_real), intent(in) :: mf       !< mass flux at the face
    integer(ccs_int), intent(in) :: bc     !< flag indicating whether cell is on boundary
    real(ccs_real), intent(out) :: coeff   !< advection coefficient to be calculated

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

  !> Calculates advection coefficient for neighbouring cell using UDS discretisation
  module subroutine calc_advection_coeff_uds(phi, mf, bc, coeff)
    type(upwind_field), intent(in) :: phi !< scalar field
    real(ccs_real), intent(in) :: mf      !< mass flux at the face
    integer(ccs_int), intent(in) :: bc    !< flag indicating whether cell is on boundary
    real(ccs_real), intent(out) :: coeff  !< advection coefficient to be calculated

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
