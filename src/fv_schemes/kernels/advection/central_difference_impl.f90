submodule(central_difference_kernel) cd_kernel_impl
  implicit none

contains
  module pure function advect_cd_eval_coeffs(self, flux_coeff) result(coeffs)
    class(cd_advection_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: flux_coeff   ! ṁ_f
    real(ccs_real), dimension(2) :: coeffs       ! (P, F)

    associate (foo => self); end associate

    coeffs(1) = max(flux_coeff, 0.0_ccs_real)   ! owner    P
    coeffs(2) = min(flux_coeff, 0.0_ccs_real)   ! neighbour F
  end function advect_cd_eval_coeffs

  module pure function advect_cd_eval_explicit(self, phi, flux_coeff, lf, rvecs, grads) result(expl)
    class(cd_advection_kernel), intent(in) :: self
    real(ccs_real), dimension(2), intent(in) :: phi
    real(ccs_real), intent(in) :: flux_coeff
    real(ccs_real), intent(in) :: lf
    real(ccs_real), dimension(3, 2), intent(in) :: rvecs
    real(ccs_real), dimension(3, 2), intent(in) :: grads
    real(ccs_real) :: expl

    real(ccs_real) :: flux_up, phi_cd
    real(ccs_real), dimension(2) :: coeffs

    associate (foo => rvecs); end associate
    associate (foo => grads); end associate

    coeffs = self%eval_coeffs(flux_coeff)
    flux_up = coeffs(1) * phi(1) + coeffs(2) * phi(2)

    ! --- central-difference face value
    phi_cd = lf * phi(1) + (1.0_ccs_real - lf) * phi(2)

    ! --- deferred correction term
    expl = flux_coeff * phi_cd - flux_up
  end function advect_cd_eval_explicit

  module pure function get_cd_width(self) result(width)
    class(cd_advection_kernel), intent(in) :: self
    integer(ccs_int) :: width
    associate (foo => self); end associate

    width = 1
  end function get_cd_width

  !> Formal order of accuracy
  module pure function get_cd_order(self) result(order)
    class(cd_advection_kernel), intent(in) :: self
    integer(ccs_int) :: order
    associate (foo => self); end associate
    order = 2
  end function get_cd_order

end submodule cd_kernel_impl
