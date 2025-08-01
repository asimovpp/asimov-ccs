submodule(central_difference_kernel) cd_kernel_impl
  implicit none

contains
!==========================================
!  Second-order central‐difference advection
!==========================================

!> Implicit coefficients for the convective flux term
!!     ṁ_f  φ_f  with  φ_f = (1 -a)*φ_P + a*φ_F
!!
!!  Returns:
!!    coeffs(1) – multiplies neighbour value  φ_F
!!    coeffs(2) – multiplies owner-cell value φ_P
!!
  module pure function advect_cd_eval_coeffs(self, flux_coeff) result(coeffs)
    class(cd_advection_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: flux_coeff   ! ṁ_f = ρ u_f A_f
    real(ccs_real), dimension(2) :: coeffs

    ! Central interpolation : φ_f = ½(φ_P + φ_F)
    coeffs(1) = self%get_interpolation_factor() * flux_coeff          ! ← neighbour-side contribution
    coeffs(2) = (1 - self%get_interpolation_factor()) * flux_coeff    ! ← owner-side contribution
  end function advect_cd_eval_coeffs

!> Explicit (deferred) part — zero for a pure CD scheme
  module pure function advect_cd_eval_explicit(self, phi, flux_coeff, lf, rvecs, grads) result(expl)
    class(cd_advection_kernel), intent(in) :: self
    real(ccs_real), dimension(2), intent(in) :: phi
    real(ccs_real), intent(in) :: flux_coeff
    real(ccs_real), intent(in) :: lf
    real(ccs_real), dimension(3, 2), intent(in) :: rvecs
    real(ccs_real), dimension(3, 2), intent(in) :: grads
    real(ccs_real) :: expl

    associate (foo => self, bar => rvecs, baz => grads, lux => lf); end associate
    associate (foo => phi); end associate
    associate (foo => flux_coeff); end associate
    expl = 0.0_ccs_real          ! nothing is deferred in a plain CD scheme
  end function advect_cd_eval_explicit

!> Stencil width (1 face on either side of the owner cell)
  module pure function get_cd_width(self) result(width)
    class(cd_advection_kernel), intent(in) :: self
    integer(ccs_int) :: width
! Silence unused-variable warnings (if you compile with -Werror)
    associate (foo => self); end associate

    width = 1
  end function get_cd_width

!> Formal order of accuracy
  module pure function get_cd_order(self) result(order)
    class(cd_advection_kernel), intent(in) :: self
    integer(ccs_int) :: order

    ! Silence unused-variable warnings (if you compile with -Werror)
    associate (foo => self); end associate

    order = 2                     ! second-order in space
  end function get_cd_order

  ! Setter
  module pure subroutine set_interpolation_factor(self, interpol_fact)
    class(cd_advection_kernel), intent(inout) :: self
    real(ccs_real), intent(in) :: interpol_fact
    self%m_interpol_factor = interpol_fact
  end subroutine set_interpolation_factor

! Getter
  module pure function get_interpolation_factor(self) result(interpol_fact)
    class(cd_advection_kernel), intent(in) :: self
    real(ccs_real) :: interpol_fact
    interpol_fact = self%m_interpol_factor
  end function get_interpolation_factor

end submodule cd_kernel_impl
