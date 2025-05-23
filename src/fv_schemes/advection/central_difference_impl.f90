submodule(fv_kernels) central_difference_advection
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
  module pure function advect_cd_eval_coeffs(self, flux_coeff, interpol_fact) result(coeffs)
    class(central_advection_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: flux_coeff   ! ṁ_f = ρ u_f A_f
    real(ccs_real), intent(in) :: interpol_fact  ! interpolation factor (0.5 for CD)
    real(ccs_real), dimension(2) :: coeffs

    ! Coefficients for the convective flux term
    if (.not. present(interpol_fact)) interpol_fact = 0.5_ccs_real ! default value

    ! Silence unused-variable warnings (if you compile with -Werror)
    associate (foo => self); end associate

    ! Central interpolation : φ_f = ½(φ_P + φ_F)
    coeffs(1) = interpol_fact * flux_coeff          ! ← neighbour-side contribution
    coeffs(2) = (1 - interpol_fact) * flux_coeff    ! ← owner-side contribution
  end function advect_cd_eval_coeffs

!> Explicit (deferred) part — zero for a pure CD scheme
  module pure function advect_cd_eval_explicit(self, flux_coeff, lf, rvecs, grads) result(expl)
    class(central_advection_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: flux_coeff
    real(ccs_real), intent(in) :: lf
    real(ccs_real), dimension(3, 2), intent(in) :: rvecs
    real(ccs_real), dimension(3, 2), intent(in) :: grads
    real(ccs_real) :: expl

    associate (foo => self, bar => rvecs, baz => grads); end associate
    expl = 0.0_ccs_real          ! nothing is deferred in a plain CD scheme
  end function advect_cd_eval_explicit

!> Stencil width (1 face on either side of the owner cell)
  module pure function get_cd_width(self) result(width)
    class(central_advection_kernel), intent(in) :: self
    integer(ccs_int) :: width
    width = 1
  end function get_cd_width

!> Formal order of accuracy
  module pure function get_cd_order(self) result(order)
    class(central_advection_kernel), intent(in) :: self
    integer(ccs_int) :: order
    order = 2                     ! second-order in space
  end function get_cd_order

end submodule advection_cds_impl
