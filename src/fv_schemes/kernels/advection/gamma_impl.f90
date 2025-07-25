submodule(fv_kernels) gamma_advection
  implicit none

contains

  module pure function advect_gamma_eval_coeffs(self, flux_coeff) result(coeffs)
    class(gamma_advection_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: flux_coeff     ! ṁ_f
    real(ccs_real), dimension(2) :: coeffs

! Silence unused-variable warnings (if you compile with -Werror)
    associate (foo => self); end associate

    coeffs(1) = max(-flux_coeff, 0.0_ccs_real)   ! neighbour F
    coeffs(2) = max(flux_coeff, 0.0_ccs_real)    ! owner    P
  end function advect_gamma_eval_coeffs

  module pure function advect_gamma_eval_explicit(self, flux_coeff, lf, rvecs, grads) result(expl)
    class(gamma_advection_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: flux_coeff        ! ρ u A
    real(ccs_real), intent(in) :: lf                ! face-interp factor (0→P,1→F)
    real(ccs_real), dimension(3, 2), intent(in) :: rvecs          ! x_Pf , x_Ff
    real(ccs_real), dimension(3, 2), intent(in) :: grads          ! ∇φ_P , ∇φ_F
    real(ccs_real) :: expl                                         ! deferred flux

    ! ---- local aliases ------------------------------------------------------
    logical :: posFlux
    real(ccs_real) :: phiP, phiF, phiUp, phiGamma
    real(ccs_real), dimension(3) :: gradP, dPF, d_up
    real(ccs_real) :: dphi, ddphi, phiPt, gamma_m
    real(ccs_real) :: w, coeffF, coeffP, bm
    real(ccs_real), dimension(2) :: phi_coeffs       ! [ φ_P , φ_F ]
    real(ccs_real), parameter :: zero = 0.0_ccs_real, one = 1.0_ccs_real

    phi_coeffs = self%eval_coeffs(flux_coeff)

    bm = self%beta_m
    posFlux = flux_coeff >= zero

    ! ---- pick UP-wind quantities depending on flow direction ---------------
    if (posFlux) then                           ! flow P → F
      phiP = phi_coeffs(1)          ! up-wind value
      phiF = phi_coeffs(2)
      gradP = grads(:, 1)
      dPF = rvecs(:, 2) - rvecs(:, 1)          ! x_F − x_P
      d_up = -rvecs(:, 1)                      ! x_f − x_P
      w = lf                               ! interpolation weight from P
    else                                        ! flow F → P
      phiP = phi_coeffs(2)
      phiF = phi_coeffs(1)
      gradP = grads(:, 2)
      dPF = rvecs(:, 1) - rvecs(:, 2)          ! x_P − x_F
      d_up = rvecs(:, 2)                      ! x_f − x_F
      w = one - lf                         ! weight from up-wind cell
    end if

    ! ---- normalised variable φ̃ --------------------------------------------
    dphi = phiF - phiP
    ddphi = 2.0_ccs_real * dot_product(gradP, dPF)

    if (abs(ddphi) < 1.0e-30_ccs_real) then
      phiPt = zero
    else
      phiPt = one - dphi / ddphi
    end if

    ! ---- choose interpolation weights according to Gamma NVD ---------------
    if (phiPt <= zero .or. phiPt >= one) then          ! → revert to UD
      coeffF = zero
      coeffP = one
    else if (phiPt > bm) then                          ! pure CDS
      coeffF = one - w
      coeffP = w
    else                                               ! bounded Gamma blend
      gamma_m = phiPt / bm
      coeffF = gamma_m * (one - w)
      coeffP = one - coeffF
    end if

    ! ---- face values and deferred correction -------------------------------
    phiGamma = coeffP * phiP + coeffF * phiF
    phiUp = phiP                                  ! 1st-order up-wind value

    expl = flux_coeff * (phiGamma - phiUp)
  end function advect_gamma_eval_explicit

  module pure function get_gamma_width(self) result(width)
    class(gamma_advection_kernel), intent(in) :: self
    integer(ccs_int) :: width

    ! Silence unused-variable warnings (if you compile with -Werror)
    associate (foo => self); end associate

    width = 1
  end function get_gamma_width

  module pure function get_gamma_order(self) result(order)
    class(gamma_advection_kernel), intent(in) :: self
    integer(ccs_int) :: order
! Silence unused-variable warnings (if you compile with -Werror)
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

end submodule gamma_advection
