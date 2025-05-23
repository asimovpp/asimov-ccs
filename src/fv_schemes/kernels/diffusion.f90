!> Implements the diffusion kernel
!  Currently there is only a single diffusion kernel implementation based on central differences.

submodule(fv_kernels) fv_diffusion_kernel

  implicit none

contains

  ! Computes the implicit coefficients for the diffusion term.
  !
  ! Note that the flux_coeff is expected of the form
  !   d_f = Gamma_f A_f / |dx|
  module pure function diffusion_coeffs(self, flux_coeff) result(coeffs)

    class(diffusion_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: flux_coeff
    real(ccs_real), dimension(2) :: coeffs

    ! Silence unused compiler warnings
    associate (foo => self)
    end associate

    coeffs(1) = flux_coeff
    coeffs(2) = -flux_coeff

  end function diffusion_coeffs

  module pure function diffusion_eval(self, flux_coeff, lf, rvecs, grads, phi_coeffs) result(expl)
    class(diffusion_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: flux_coeff
    real(ccs_real), intent(in) :: lf
    real(ccs_real), dimension(3, 2), intent(in) :: rvecs
    real(ccs_real), dimension(3, 2), intent(in) :: grads
    real(ccs_real), dimension(2), intent(in), optional :: phi_coeffs
    real(ccs_real) :: expl

    ! Silence unused compiler warnings
    associate (foo => self)
    end associate
    associate (foo => lf, bar => rvecs, baz => grads, qux => phi_coeffs)
    end associate

    expl = 0.0_ccs_real
    expl = flux_coeff * expl

  end function diffusion_eval

  module pure function diffusion_width(self) result(width)
    class(diffusion_kernel), intent(in) :: self
    integer(ccs_int) :: width

    ! Silence unused compiler warnings
    associate (foo => self)
    end associate

    width = 1

  end function diffusion_width

  module pure function diffusion_order(self) result(order)
    class(diffusion_kernel), intent(in) :: self
    integer(ccs_int) :: order

    ! Silence unused compiler warnings
    associate (foo => self)
    end associate

    order = 2

  end function diffusion_order

end submodule fv_diffusion_kernel
