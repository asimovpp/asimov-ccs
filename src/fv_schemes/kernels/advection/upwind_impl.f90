submodule(upwind_kernel) upwind_impl
  implicit none

contains
!> Calculates advection coefficient for neighbouring cell using UDS discretisation
  module pure function advect_upwind_eval_coeffs(self, flux_coeff) result(coeffs)
    class(upwind_advection_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: flux_coeff
    real(ccs_real), dimension(2) :: coeffs       ! (P, F)

    associate (foo => self); end associate

    coeffs(1) = max(flux_coeff, 0.0_ccs_real)   ! owner    P
    coeffs(2) = min(flux_coeff, 0.0_ccs_real)   ! neighbour F
  end function advect_upwind_eval_coeffs

  module pure function advect_upwind_eval_explicit(self, phi, flux_coeff, lf, rvecs, grads) result(expl)
    class(upwind_advection_kernel), intent(in) :: self
    real(ccs_real), dimension(2), intent(in) :: phi
    real(ccs_real), intent(in) :: flux_coeff
    real(ccs_real), intent(in) :: lf
    real(ccs_real), dimension(3, 2), intent(in) :: rvecs
    real(ccs_real), dimension(3, 2), intent(in) :: grads

    real(ccs_real) :: expl

    ! Silence unused compiler warnings
    associate (foo => self, bar => rvecs, baz => grads, lux => lf, psi => flux_coeff); end associate
    associate (foo => phi); end associate

    ! Calculate explicit term
    expl = 0.0_ccs_real

  end function advect_upwind_eval_explicit

  module pure function get_upwind_width(self) result(width)
    class(upwind_advection_kernel), intent(in) :: self
    integer(ccs_int) :: width

! Silence unused-variable warnings (if you compile with -Werror)
    associate (foo => self); end associate

    width = 1
  end function get_upwind_width

  module pure function get_upwind_order(self) result(order)
    class(upwind_advection_kernel), intent(in) :: self
    integer(ccs_int) :: order

    ! Silence unused-variable warnings (if you compile with -Werror)
    associate (foo => self); end associate

    order = 1
  end function get_upwind_order

end submodule upwind_impl
