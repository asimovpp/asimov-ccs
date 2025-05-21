program test_diffusion_kernel
  !! Test the diffusion kernel by running a refinement loop, computing the discretisation
  !!  error (difference between analytical and computed values), and checking the order of
  !!  convergence is the theoretical order for the scheme

  use testing_lib 

  use kinds, only: ccs_real, ccs_int
  use error_analysis, only: compute_order

  use fv_kernels, only: diffusion_kernel

  implicit none

  integer, parameter :: num_iters = 10
  integer(ccs_int) :: i
  real(ccs_real) :: phi_P, interpol_factor
  real(ccs_real) :: phi_N
  real(ccs_real) :: dx
  real(ccs_real), dimension(3) :: x_N, x_P, x_f
  real(ccs_real), dimension(2) :: coeffs
  real(ccs_real) :: rhs
  real(ccs_real) :: mu_f
  real(ccs_real) :: flux_coeff
  real(ccs_real), dimension(num_iters) :: errors
  real(ccs_real), dimension(num_iters) :: refinements
  real(ccs_real) :: order

  real(ccs_real), dimension(3, 2) :: rvecs, grads

  type(diffusion_kernel) :: diffusion

  ! Initialising values
  x_f = [0, 0, 0]
  ! interpol_factor = 2.0_ccs_real / 3.0_ccs_real
  interpol_factor = 0.5_ccs_real
  rvecs = 0.0_ccs_real
  grads = 0.0_ccs_real

  ! Refinement loop
  do i = 1, num_iters
    dx = 1.0e-1_ccs_real / 2**i
    x_P = x_f - [(1.0_ccs_real - interpol_factor) * dx, 0.0_ccs_real, 0.0_ccs_real]
    x_N = x_f + [interpol_factor * dx, 0.0_ccs_real, 0.0_ccs_real]

    phi_P = phi(x_P(1))
    phi_N = phi(dx)
    mu_f = interpol_factor * mu(x_P(1)) + (1 - interpol_factor) * mu(x_N(1))

    flux_coeff = mu_f / dx
    coeffs = diffusion%eval_coeffs(flux_coeff)
    rhs = diffusion%eval_explicit(flux_coeff, interpol_factor, rvecs, grads)

    ! Compute discretisation error
    call get_error(coeffs, rhs, x_P, x_N, x_f, errors(i))
    refinements(i) = dx
  end do

  ! Compute convergence order
  call compute_order(refinements, errors, order)

  ! XXX: Diffusion is only 2nd order accurate for equally spaced grids (i.e. interpol_factor=0.5)
  call assert_gt(order, diffusion%get_order() * 0.95_ccs_real, "Convergence order not preserved by diffusion kernel")

contains

  real(ccs_real) pure function mu(x)
    ! Velocity field
    real(ccs_real), intent(in) :: x

    ! mu = cos(x)
    mu = 3.1415_ccs_real + 0 * x
  end function mu

  real(ccs_real) pure function phi(x)
    ! Function being tested
    real(ccs_real), intent(in) :: x

    phi = sin(2 * x + 17)
    ! phi = 2.0_ccs_real * x
  end function phi

  real(ccs_real) pure function dphi(x)
    ! Derivative of function being tested
    real(ccs_real), intent(in) :: x

    dphi = 2 * cos(2 * x + 17)
    ! dphi = 2.0_ccs_real + 0 * x
  end function dphi

  subroutine get_error(coeffs, rhs, x_P, x_N, x_f, error)
    ! Computes discretisation error by comparing with analytical value 
    real(ccs_real), intent(in) :: coeffs(2), rhs
    real(ccs_real), intent(in) :: x_P(3)
    real(ccs_real), intent(in) :: x_N(3)
    real(ccs_real), intent(in) :: x_f(3)
    real(ccs_real), intent(out) :: error
    real(ccs_real) :: analytical

    analytical = mu_f * dphi(x_f(1))
    print *, analytical, coeffs(1) * phi(x_P(1)), coeffs(2) * phi(x_N(1)), rhs, rhs - (coeffs(1) * phi(x_P(1)) + coeffs(2) * phi(x_N(1)))
    error = abs(analytical - (rhs - (coeffs(1) * phi(x_P(1)) + coeffs(2) * phi(x_N(1)))))
  end subroutine get_error


end program test_diffusion_kernel
