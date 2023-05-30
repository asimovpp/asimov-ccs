!v Submodule file fv_discretisation.smod
!
!  Implementations of the finite volume method using the various discretisation schemes scheme
!
!  @build discretisation

submodule(fv) fv_discretisation

  use meshing, only: get_face_interpolation
  implicit none

contains

  !> Calculates advection coefficient for neighbouring cell using CDS discretisation
  module subroutine calc_advection_coeff_cds(phi, loc_f, mf, bc, coeff)
    type(central_field), intent(in) :: phi !< scalar field
    type(face_locator), intent(in) :: loc_f !< face locator
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
      call get_face_interpolation(loc_f, interpolation_factor)
      interpolation_factor = 1.0_ccs_real - interpolation_factor
    else
      interpolation_factor = 0.5_ccs_real !1.0_ccs_real
    end if
    coeff = interpolation_factor
  end subroutine calc_advection_coeff_cds

  !> Calculates advection coefficient for neighbouring cell using UDS discretisation
  module subroutine calc_advection_coeff_uds(phi, loc_f, mf, bc, coeff)
    type(upwind_field), intent(in) :: phi !< scalar field
    type(face_locator), intent(in) :: loc_f !< face locator
    real(ccs_real), intent(in) :: mf      !< mass flux at the face
    integer(ccs_int), intent(in) :: bc    !< flag indicating whether cell is on boundary
    real(ccs_real), intent(out) :: coeff  !< advection coefficient to be calculated

    ! Dummy usage to prevent unused argument.
    associate (scalar => phi, foo => bc, bar => loc_f)
    end associate

    if (mf < 0.0) then
      coeff = 1.0_ccs_real
    else
      coeff = 0.0_ccs_real
    end if

  end subroutine calc_advection_coeff_uds

end submodule fv_discretisation
