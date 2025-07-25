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

  ! Computes the explicit contribution of the diffusion term,
  module pure function diffusion_eval(self, flux_coeff, lf, rvecs, grads) result(expl)

    class(diffusion_kernel), intent(in) :: self          !< The diffusion kernel object
    real(ccs_real), intent(in) :: flux_coeff             !< The flux coefficient: Gamma / dx * A
    real(ccs_real), intent(in) :: lf                     !< Cell-face linear interpolatino coefficient
    real(ccs_real), dimension(3, 2), intent(in) :: rvecs !< Cell-face vectors [r_{Pf}, r_{Ff}]
    real(ccs_real), dimension(3, 2), intent(in) :: grads !< Solution gradients
    real(ccs_real):: expl                                !< The resulting explicit contribution

    ! Silence unused compiler warnings
    associate(foo => self, bar => lf)
    end associate

    ! Non-orthogonality correction (diffusive flux) (Ferziger & Peric 4th ed, sec 9.7.2)
    expl = (dot_product(grads(:, 2), rvecs(:, 2)) - dot_product(grads(:, 1), rvecs(:, 1)))
    expl = -flux_coeff * expl
    
  end function diffusion_eval

  ! Returns the diffusion stencil width.
  module pure function diffusion_width(self) result(width)
    class(diffusion_kernel), intent(in) :: self
    integer(ccs_int) :: width

    ! Silence unused compiler warnings
    associate (foo => self)
    end associate

    width = 1

  end function diffusion_width

  ! Returms the expected order of the discretisation.
  module pure function diffusion_order(self) result(order)
    class(diffusion_kernel), intent(in) :: self
    integer(ccs_int) :: order

    ! Silence unused compiler warnings
    associate (foo => self)
    end associate

    order = 2

  end function diffusion_order

end submodule fv_diffusion_kernel
