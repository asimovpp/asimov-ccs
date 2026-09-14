submodule(linear_upwind_kernel) luds_impl
  implicit none

contains

  module pure function advect_luds_eval_coeffs(self, flux_coeff) result(coeffs)
    class(luds_advection_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: flux_coeff   ! rho u A with sign
    real(ccs_real), dimension(2) :: coeffs       ! (P, F)

    associate (foo => self); end associate

    coeffs(1) = max(flux_coeff, 0.0_ccs_real)   ! owner    P
    coeffs(2) = min(flux_coeff, 0.0_ccs_real)   ! neighbour F
  end function advect_luds_eval_coeffs

  module pure function advect_luds_eval_explicit(self, phi, flux_coeff, lf, rvecs, grads) result(expl)
    class(luds_advection_kernel), intent(in) :: self
    real(ccs_real), dimension(2), intent(in) :: phi
    real(ccs_real), intent(in) :: flux_coeff      ! rho u A
    real(ccs_real), intent(in) :: lf              ! unused here
    real(ccs_real), dimension(3, 2), intent(in) :: rvecs        ! r_Pf , r_Ff
    real(ccs_real), dimension(3, 2), intent(in) :: grads        ! grad phi_P, grad phi_F
    real(ccs_real) :: expl                                     ! high-order - UD

    ! ---------------- bookkeeping ------------------------------------
    logical :: pos
    real(ccs_real) :: phi_lud
    real(ccs_real), dimension(3) :: grad_up, d_up

    associate (foo => self); end associate
    associate (foo => phi); end associate
    associate (foo => lf); end associate

    pos = flux_coeff >= 0.0_ccs_real
    if (pos) then                       ! flow P -> F
      grad_up = grads(:, 1)
      d_up = rvecs(:, 1)               ! x_f - x_P
    else                                ! flow F -> P
      grad_up = grads(:, 2)
      d_up = rvecs(:, 2)                ! x_f - x_F
    end if

    phi_lud = dot_product(grad_up, d_up)   ! linear reconstruction

    ! LUDS computes face value as upwind + linear extrapolation, therefore
    ! the deferred correction (f_LUDS - f_U) reduces to the linear extrapolation.
    expl = flux_coeff * phi_lud         ! deferred correction
  end function advect_luds_eval_explicit

  module pure function get_luds_width(self) result(width)
    class(luds_advection_kernel), intent(in) :: self
    integer(ccs_int) :: width
    associate (foo => self); end associate

    width = 1
  end function get_luds_width

  module pure function get_luds_order(self) result(order)
    class(luds_advection_kernel), intent(in) :: self
    integer(ccs_int) :: order

    associate (foo => self); end associate

    order = 2
  end function get_luds_order

end submodule luds_impl
