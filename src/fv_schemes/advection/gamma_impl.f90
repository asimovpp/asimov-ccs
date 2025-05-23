submodule(fv_kernels) gamma_advection
  implicit none

contains

  module pure function advect_gamma_eval_coeffs(self, flux_coeff) result(coeffs)
    class(gamma_advection_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: flux_coeff     ! ṁ_f
    real(ccs_real), dimension(2) :: coeffs

    coeffs(1) = max(-flux_coeff, 0.0_ccs_real)   ! neighbour F
    coeffs(2) = max(flux_coeff, 0.0_ccs_real)    ! owner    P
  end function advect_gamma_eval_coeffs

  module pure function advect_gamma_eval_explicit(self, flux_coeff, lf, rvecs, grads) result(expl)
    class(gamma_advection_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: flux_coeff           ! ṁ_f
    real(ccs_real), intent(in) :: lf                   ! f_x interpolation factor
    real(ccs_real), dimension(3, 2), intent(in) :: rvecs             ! vectors P→f and F→f
    real(ccs_real), dimension(3, 2), intent(in) :: grads             ! ∇φ_P , ∇φ_F
    real(ccs_real) :: expl                                          ! deferred flux

    !--- local aliases ---------------------------------------------------------
    real(ccs_real) :: phiP, phiF            ! scalar values in owner / neighbour
    real(ccs_real) :: phiPt, gamma_m, coeffF, coeffP, phiGamma, phiUpwind
    real(ccs_real), dimension(3) :: d
    real(ccs_real) :: dphi, ddphi
    real(ccs_real), parameter :: one = 1.0_ccs_real
    real(ccs_real) :: bm
    bm = self % beta_m
    ! NOTE: the caller must guarantee that phi values live at the same indices
    !       as the gradients just passed in, i.e.  grads(:,1) => P, grads(:,2) => F
    phiP = lf        ! <-- convention: caller passes φ_P  in lf   (see below)
    phiF = rvecs(1, 1)   ! <-- convention: caller packs φ_F in rvecs(1,1)

    !--- build the NVD normalised variable  φ~ = φPt ---------------------------
    d = rvecs(:, 1) - rvecs(:, 2)             ! d = xF - xP   (owner→neigh.)
    dphi = phiF - phiP
    ddphi = 2.0_ccs_real * dot_product(grads(:, 1), d)   ! 2 ∇φ_P · d
    phiPt = one - dphi / max(ddphi, 1.0e-30_ccs_real)   ! avoid /0

    !--- choose weights for Gamma scheme --------------------------------------
    if (phiPt <= 0.0_ccs_real .or. phiPt >= 1.0_ccs_real) then
      ! Outside (0,1)  →  revert to first-order up-wind
      coeffF = 0.0_ccs_real
      coeffP = 1.0_ccs_real

    else if (phiPt > self % beta_m) then
      ! Inside (β_m , 1)  →  pure central difference
      coeffF = 1.0_ccs_real - lf
      coeffP = lf

    else
      ! 0 < φ̃ ≤ β_m  →  Gamma blend
      gamma_m = phiPt / self % beta_m
      coeffF = gamma_m * (1.0_ccs_real - lf)
      coeffP = 1.0_ccs_real - coeffF
    end if

    !--- high-order flux minus up-wind flux -----------------------------------
    phiGamma = coeffP * phiP + coeffF * phiF
    if (flux_coeff >= 0.0_ccs_real) then
      phiUpwind = phiP
    else
      phiUpwind = phiF
    end if

    expl = flux_coeff * (phiGamma - phiUpwind)   ! deferred correction
  end function advect_gamma_eval_explicit

  module pure function get_gamma_width(self) result(width)
    class(gamma_advection_kernel), intent(in) :: self
    integer(ccs_int) :: width
    width = 1
  end function get_gamma_width

  module pure function get_gamma_order(self) result(order)
    class(gamma_advection_kernel), intent(in) :: self
    integer(ccs_int) :: order
    order = 2
  end function get_gamma_order

  module pure subroutine set_beta_m(self, new_bm)
    class(gamma_advection_kernel), intent(inout) :: self
    real(ccs_real), intent(in) :: new_bm
    self % beta_m = max(0.10_ccs_real, min(0.50_ccs_real, new_bm))
  end subroutine set_beta_m

  module pure function get_beta_m(self) result(bm)
    class(gamma_advection_kernel), intent(in) :: self
    real(ccs_real) :: bm
    bm = self % beta_m
  end function get_beta_m

end submodule gamma_advection
