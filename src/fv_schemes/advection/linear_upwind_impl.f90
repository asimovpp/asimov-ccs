submodule(fv_kernels) linear_upwind_advection
  implicit none

contains

  module pure function advect_luds_eval_coeffs(self, flux_coeff) result(coeffs)
    class(luds_advection_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: flux_coeff   ! ρ u A with sign
    real(ccs_real), dimension(2) :: coeffs       ! (F, P)

    coeffs(1) = max(-flux_coeff, 0.0_ccs_real)   ! neighbour
    coeffs(2) = max(flux_coeff, 0.0_ccs_real)   ! owner
  end function advect_luds_eval_coeffs

  module pure function advect_luds_eval_explicit(self, flux_coeff, lf, rvecs, grads, phi_coeffs) result(expl)
    class(luds_advection_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: flux_coeff      ! ρ u A
    real(ccs_real), intent(in) :: lf              ! unused here
    real(ccs_real), dimension(3, 2), intent(in) :: rvecs        ! r_Pf , r_Ff
    real(ccs_real), dimension(3, 2), intent(in) :: grads        ! ∇φ_P , ∇φ_F
    real(ccs_real), dimension(2), intent(in), optional :: phi_coeffs     ! [ φ_P , φ_F ]
    real(ccs_real) :: expl                                     ! high-order – UD

  !! ---------- up-wind / down-wind bookkeeping ---------------------------
    logical :: posFlux
    real(ccs_real) :: phiP, phiF
    real(ccs_real), dimension(3) :: gradP, d_pf, d_up
    real(ccs_real) :: dphi, ddphi, phiPt, phiLUDS, phiUp
    real(ccs_real), parameter :: one = 1.0_ccs_real

    posFlux = flux_coeff >= 0.0_ccs_real

    if (posFlux) then
      ! flux P → F  :  owner  is up-wind
      phiP = phi_coeffs(1); gradP = grads(:, 1)
      phiF = phi_coeffs(2)
      d_pf = rvecs(:, 1) - rvecs(:, 2)    ! x_F − x_P
      d_up = rvecs(:, 1)                 ! x_f − x_P
    else
      ! flux F → P  :  neighbour is up-wind
      phiP = phi_coeffs(2); gradP = grads(:, 2)
      phiF = phi_coeffs(1)
      d_pf = -(rvecs(:, 1) - rvecs(:, 2)) ! x_P − x_F
      d_up = -rvecs(:, 2)                ! x_f − x_F
    end if

  !! ---------- boundedness check (Ferziger & Perić, §5.5) ---------------
    dphi = phiF - phiP
    ddphi = 2.0_ccs_real * dot_product(gradP, d_pf)

    if (ddphi == 0.0_ccs_real) then
      phiPt = 0.0_ccs_real          ! degenerate → fall back to UD
    else
      phiPt = one - dphi / ddphi
    end if

    if (phiPt <= 0.0_ccs_real .or. phiPt >= 1.0_ccs_real) then
      ! ---------- out of bounds → pure up-wind ----------------------------
      phiLUDS = phiP
    else
      ! ---------- genuine second-order linear up-wind ---------------------
      phiLUDS = phiP + dot_product(gradP, d_up)
    end if

    phiUp = phiP                    ! 1st-order up-wind value
    expl = flux_coeff * (phiLUDS - phiUp)
  end function advect_luds_eval_explicit

  module pure function get_luds_width(self) result(width)
    class(luds_advection_kernel), intent(in) :: self
    integer(ccs_int) :: width
    width = 1                      ! compact stencil
  end function get_luds_width

  module pure function get_luds_order(self) result(order)
    class(luds_advection_kernel), intent(in) :: self
    integer(ccs_int) :: order
    order = 2                      ! formal second order
  end function get_luds_order

end submodule linear_upwind_advection
