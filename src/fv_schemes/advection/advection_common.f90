submodule(fv_kernels) advection_common
  use kinds, only: ccs_real
  implicit none

contains

  !> Base implementation of the advection kernel
  !> Simple prototype
  module pure function advection_eval_coeffs(self, flux_coeff) result(coeffs)
    class(advection_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: flux_coeff
    real(ccs_real), dimension(2) :: coeffs

    ! Silence unused dummy arguments
    associate(foo=>self); end associate
    associate(foo=>flux_coeff); end associate

    coeffs = [1.0_ccs_real, 0.0_ccs_real]  ! Example coefficients
  end function advection_eval_coeffs

  module pure function advection_eval_explicit(self, flux_coeff, lf, rvecs, grads, phi_coeffs) result(expl)
    class(advection_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: flux_coeff
    real(ccs_real), intent(in) :: lf
    real(ccs_real), dimension(3, 2), intent(in) :: rvecs
    real(ccs_real), dimension(3, 2), intent(in) :: grads
    real(ccs_real), dimension(2), intent(in), optional :: phi_coeffs
    real(ccs_real) :: expl

    ! Silence unused dummy arguments
    associate(foo=>self); end associate
    associate(foo=>flux_coeff); end associate
    associate(foo=>lf); end associate
    associate(foo=>rvecs); end associate
    associate(foo=>grads); end associate
    associate(foo=>phi_coeffs); end associate
      
    expl = 0.0_ccs_real  ! Example explicit term
  end function advection_eval_explicit

  module pure function advection_width(self) result(width)
    class(advection_kernel), intent(in) :: self
    integer(ccs_int) :: width

    ! Silence unused dummy arguments
    associate(foo=>self); end associate
      
    width = 1
  end function advection_width

  module pure function advection_order(self) result(order)
    class(advection_kernel), intent(in) :: self
    integer(ccs_int) :: order

    ! Silence unused dummy arguments
    associate(foo=>self); end associate

    order = 1
  end function advection_order

end submodule advection_common
