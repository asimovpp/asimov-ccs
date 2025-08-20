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
    real(ccs_real), intent(in) :: lf                ! face-interp factor (1→P,0→F)
    real(ccs_real), dimension(3, 2), intent(in) :: rvecs          ! x_Pf , x_Ff
    real(ccs_real), dimension(3, 2), intent(in) :: grads          ! ∇φ_P , ∇φ_F
    real(ccs_real) :: expl                                         ! deferred flux

    ! ---------------- parameters & helpers ---------------------------
    real(ccs_real) :: bm, phi_up, phi_dn, phi_pt, phi_g
    real(ccs_real) :: phi_cd, gamma_m
    real(ccs_real), dimension(3) :: grad_up, d_PF
    !
    bm = self%beta_m
    call get_flow_orientation(flux_coeff, phi, grads, rvecs, phi_up, phi_dn, grad_up, d_PF)

    ! Central difference
    phi_cd = lf * phi(1) + (1.0_ccs_real - lf) * phi(2)

    ! normalised variable
    phi_pt = compute_nvd(phi_up, phi_dn, grad_up, d_PF)

    ! Gamma NVD
    gamma_m = compute_gamma(phi_pt, bm)
    phi_g = gamma_m * phi_cd + (1.0_ccs_real - gamma_m) * phi_up

    expl = flux_coeff * (phi_g - phi_up)
  end function advect_gamma_eval_explicit

  pure subroutine get_flow_orientation(flux_coeff, phi, grads, rvecs, phi_C, phi_D, grad_C, d)
    real(ccs_real), intent(in) :: flux_coeff
    real(ccs_real), dimension(2), intent(in) :: phi
    real(ccs_real), dimension(3, 2), intent(in) :: grads
    real(ccs_real), dimension(3, 2), intent(in) :: rvecs
    real(ccs_real), intent(out) :: phi_C
    real(ccs_real), intent(out) :: phi_D
    real(ccs_real), dimension(3), intent(out) :: grad_C
    real(ccs_real), dimension(3), intent(out) :: d

    if (flux_coeff >= 0.0_ccs_real) then
      phi_C = phi(1)
      phi_D = phi(2)
      grad_C = grads(:, 1)
      d = rvecs(:, 1) - rvecs(:, 2)
    else                       
      phi_C = phi(2)
      phi_D = phi(1)
      grad_C = grads(:, 2)
      d = rvecs(:, 2) - rvecs(:, 1) 
    end if
  end subroutine get_flow_orientation

  real(ccs_real) pure function compute_nvd(phi_C, phi_D, grad_C, d)
    real(ccs_real), intent(in) :: phi_C, phi_D
    real(ccs_real), dimension(3), intent(in) :: grad_C
    real(ccs_real), dimension(3), intent(in) :: d

    real(ccs_real) :: d_grad_C

    d_grad_C = dot_product(grad_C, d)

    if (abs(d_grad_C) < 100.0_ccs_real * tiny(1.0_ccs_real)) then
      compute_nvd = 0.0_ccs_real
    else
      compute_nvd = 1.0_ccs_real - (phi_D - phi_C) / (2.0_ccs_real * d_grad_C)
    end if
    
  end function compute_nvd

  real(ccs_real) pure function compute_gamma(phi_nvd, beta)
    real(ccs_real), intent(in) :: phi_nvd
    real(ccs_real), intent(in) :: beta
    
    if (phi_nvd <= 0.0_ccs_real .or. phi_nvd >= 1.0_ccs_real) then
      ! Upwind
      compute_gamma = 0.0_ccs_real
    else if (phi_nvd > beta) then
      ! Central
      compute_gamma = 1.0_ccs_real
    else
      ! Blended
      compute_gamma = phi_nvd / beta
    end if

  end function compute_gamma

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
    self%beta_m = max(0.10_ccs_real, min(0.50_ccs_real, new_bm))
  end subroutine set_beta_m

  module pure function get_beta_m(self) result(bm)
    class(gamma_advection_kernel), intent(in) :: self
    real(ccs_real) :: bm
    bm = self%beta_m
  end function get_beta_m

end submodule gamma_impl
