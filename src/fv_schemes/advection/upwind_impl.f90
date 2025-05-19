submodule(fv_kernels:advection_common) upwind_advection
implicit none

contains
!> Calculates advection coefficient for neighbouring cell using UDS discretisation
module pure function advect_upwind_coeffs(self, flux_coeff) result(coeffs)
  class(upwind_advection_kernel), intent(in) :: self
  real(ccs_real), intent(in) :: flux_coeff
  real(ccs_real), dimension(2) :: coeffs
  real(ccs_real) :: a_P  !< advection coefficient for current cell
  real(ccs_real) :: a_F  !< advection coefficient for neighbour cell

  ! Silence unused compiler warnings
  associate (foo => self)
  end associate

  ! Calculate coefficients
  a_P  = max( F_f, 0.0_ccs_real )
  a_F  = max(-F_f, 0.0_ccs_real )
  coeffs = [a_P, a_F]
end function advect_upwind_coeffs

module subroutine advect_upwind_eval(self, flux_coeff, lf, rvecs, grads) result(expl)
  class(upwind_advection_kernel), intent(in) :: self
  real(ccs_real), intent(in) :: flux_coeff
  real(ccs_real), intent(in) :: lf
  real(ccs_real), dimension(3, 2), intent(in) :: rvecs
  real(ccs_real), dimension(3, 2), intent(in) :: grads
  real(ccs_real) :: expl

  ! Silence unused compiler warnings
    associate(foo => self)
    end associate
    associate(foo => lf, bar => rvecs, baz => grads)
    end associate

  ! Calculate explicit term
  expl = 0.0_ccs_real

end subroutine advect_upwind_eval

module pure function get_upwind_width(self) result(width)
  class(upwind_advection_kernel), intent(in) :: self
  integer(ccs_int) :: width
  width = 1
end function get_upwind_width

module pure function get_upwind_order(self) result(order)
  class(upwind_advection_kernel), intent(in) :: self
  integer(ccs_int) :: order
  order = 1
end function get_upwind_order

end submodule advection_upwind
