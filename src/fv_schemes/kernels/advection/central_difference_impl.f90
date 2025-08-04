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
    real(ccs_real), intent(in)              :: flux_coeff   ! ṁ_f
    real(ccs_real), dimension(2)            :: coeffs       ! (F , P)

    associate (foo => self); end associate

    coeffs(1) = min(flux_coeff, 0.0_ccs_real)   ! neighbour F
    coeffs(2) = max(flux_coeff, 0.0_ccs_real)   ! owner    P
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

    real(ccs_real), dimension(2) :: coeffs
    real(ccs_real)               :: wP, wF, lambda_eff
    real(ccs_real)               :: phi_up, phi_dn, phi_cds

    associate (foo => rvecs); end associate
    associate (foo => grads); end associate

    coeffs = self%eval_coeffs(flux_coeff)
    wP     = 0.5_ccs_real * ( 1.0_ccs_real + sign(1.0_ccs_real, flux_coeff))      ! 1 for +flow, 0 for –flow
    wF     = 1.0_ccs_real - wP                    ! complementary

    phi_up = wP*phi(1) + wF*phi(2)
    phi_dn = wF*phi(1) + wP*phi(2)                ! opposite cell

    ! --- effective interpolation factor -----------------------------
    lambda_eff = wP*lf + wF*(1.0_ccs_real - lf)   ! λ  or 1–λ

    ! --- CDS face value and deferred term ---------------------------
    phi_cds = phi_up + lambda_eff*(phi_dn - phi_up)
    expl    = flux_coeff * (phi_cds - phi_up)
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
