submodule(gamma_kernel) gamma_impl
  implicit none

contains

  module pure function advect_gamma_eval_coeffs(self, flux_coeff) result(coeffs)
    class(gamma_advection_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: flux_coeff     ! ṁ_f
    real(ccs_real), dimension(2) :: coeffs       ! (P, F)

    associate (foo => self); end associate

    coeffs(1) = max(flux_coeff, 0.0_ccs_real)   ! owner    P
    coeffs(2) = min(flux_coeff, 0.0_ccs_real)   ! neighbour F
  end function advect_gamma_eval_coeffs

  module pure function advect_gamma_eval_explicit(self, phi, flux_coeff, lf, rvecs, grads) result(expl)
    class(gamma_advection_kernel), intent(in) :: self
    real(ccs_real), dimension(2), intent(in) :: phi ! φ_F , φ_P
    real(ccs_real), intent(in) :: flux_coeff        ! ρ u A
    real(ccs_real), intent(in) :: lf                ! face-interp factor (0→P,1→F)
    real(ccs_real), dimension(3, 2), intent(in) :: rvecs          ! x_Pf , x_Ff
    real(ccs_real), dimension(3, 2), intent(in) :: grads          ! ∇φ_P , ∇φ_F
    real(ccs_real) :: expl                                         ! deferred flux

    ! ---------------- parameters & helpers ---------------------------
    real(ccs_real) :: bm, w, phi_up, phi_dn, phi_pt
    real(ccs_real) :: phi_CDS, gamma_m
    real(ccs_real), dimension(3) :: grad_up, d_PF, d_up
    logical :: pos
    !
    bm = self % beta_m
    pos = flux_coeff >= 0.0_ccs_real
    if (pos) then                 ! P → F
      phi_up = phi(1); phi_dn = phi(2)
      grad_up = grads(:, 1)
      d_PF = rvecs(:, 2) - rvecs(:, 1)
      d_up = -rvecs(:, 1)
      w = lf
    else                          ! F → P
      phi_up = phi(2); phi_dn = phi(1)
      grad_up = grads(:, 2)
      d_PF = rvecs(:, 1) - rvecs(:, 2)
      d_up = rvecs(:, 2)
      w = 1.0_ccs_real - lf
    end if
    ! normalised variable
    if (abs(dot_product(grad_up, d_PF)) < 100.0_ccs_real * tiny(1.0_ccs_real)) then
      phi_pt = 0.0_ccs_real
    else
      phi_pt = 1.0_ccs_real - (phi_dn - phi_up) / (2.0_ccs_real * dot_product(grad_up, d_PF))
    end if
    ! Gamma NVD
    if (phi_pt <= 0.0_ccs_real .or. phi_pt >= 1.0_ccs_real) then          ! → UD
      phi_CDS = phi_up
    else if (phi_pt > bm) then                          ! pure CDS
      phi_CDS = (1.0_ccs_real - w) * phi_dn + w * phi_up
    else                                               ! blended
      gamma_m = phi_pt / bm
      phi_CDS = gamma_m * ((1.0_ccs_real - w) * phi_dn) + (1.0_ccs_real - gamma_m * w) * phi_up
    end if
    expl = flux_coeff * (phi_CDS - phi_up)
  end function advect_gamma_eval_explicit

  module pure function get_gamma_width(self) result(width)
    class(gamma_advection_kernel), intent(in) :: self
    integer(ccs_int) :: width

    associate (foo => self); end associate

    width = 1
  end function get_gamma_width

  module pure function get_gamma_order(self) result(order)
    class(gamma_advection_kernel), intent(in) :: self
    integer(ccs_int) :: order
    associate (foo => self); end associate

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

end submodule gamma_impl
