submodule(fv_kernels) linear_upwind_advection
  implicit none

contains

  module pure function advect_luds_eval_coeffs(self, flux_coeff) result(coeffs)
    class(luds_advection_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: flux_coeff   ! ρ u A with sign
    real(ccs_real), dimension(2) :: coeffs       ! (F, P)

    ! Silence unused-variable warnings (if you compile with -Werror)
    associate (foo => self); end associate

    coeffs(1) = min( flux_coeff, 0.0_ccs_real)   ! neighbour F  (negative if flow P←F)
    coeffs(2) = max( flux_coeff, 0.0_ccs_real)   ! owner    P  (positive if flow P→F)
  end function advect_luds_eval_coeffs

  module pure function advect_luds_eval_explicit(self, phi, flux_coeff, lf, rvecs, grads) result(expl)
    class(luds_advection_kernel), intent(in) :: self
    real(ccs_real), dimension(2), intent(in) :: phi
    real(ccs_real), intent(in) :: flux_coeff      ! ρ u A
    real(ccs_real), intent(in) :: lf              ! unused here
    real(ccs_real), dimension(3, 2), intent(in) :: rvecs        ! r_Pf , r_Ff
    real(ccs_real), dimension(3, 2), intent(in) :: grads        ! ∇φ_P , ∇φ_F
    real(ccs_real) :: expl                                     ! high-order – UD
  ! ---------------- bookkeeping ------------------------------------
  logical                      :: pos
  real(ccs_real)               :: phiUp, phiHO
  real(ccs_real), dimension(3) :: gradUp, d_up

  associate (foo => self); end associate
  associate (foo => lf); end associate

  pos = flux_coeff >= 0.0_ccs_real
  if (pos) then                      ! flow P → F
     phiUp  = phi(1)
     gradUp = grads(:,1)
     d_up   = -rvecs(:,1)            ! x_f − x_P
  else                                ! flow F → P
     phiUp  = phi(2)
     gradUp = grads(:,2)
     d_up   =  rvecs(:,2)            ! x_f − x_F
  end if

  phiHO = phiUp + dot_product(gradUp, d_up)   ! linear reconstruction
  expl  = flux_coeff * (phiHO - phiUp)        ! deferred correction
  end function advect_luds_eval_explicit

  module pure function get_luds_width(self) result(width)
    class(luds_advection_kernel), intent(in) :: self
    integer(ccs_int) :: width
! Silence unused-variable warnings (if you compile with -Werror)
    associate (foo => self); end associate

    width = 1
  end function get_luds_width

  module pure function get_luds_order(self) result(order)
    class(luds_advection_kernel), intent(in) :: self
    integer(ccs_int) :: order

! Silence unused-variable warnings (if you compile with -Werror)
    associate (foo => self); end associate

    order = 2                      ! formal second order
  end function get_luds_order

end submodule linear_upwind_advection
