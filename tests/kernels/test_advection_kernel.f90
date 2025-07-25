program test_advection_kernel
    !! Test the advection kernel by running a refinement loop, computing the discretisation
    !!  error (difference between analytical and computed values), and checking the order of
    !!  convergence is the theoretical order for the scheme

    use testing_lib 

    use kinds, only: ccs_real, ccs_int
    use error_analysis, only: compute_order
    use fv_kernels

    implicit none

    integer(ccs_int), parameter :: num_iters = 10
    
    integer(ccs_int) :: i
    real(ccs_real) :: phi_P, interpol_factor
    real(ccs_real) :: phi_N
    real(ccs_real) :: dx
    real(ccs_real), dimension(3) :: x_N, x_P, x_f
    real(ccs_real), dimension(2) :: coeffs
    real(ccs_real) :: rhs
    real(ccs_real) :: u_f
    real(ccs_real), dimension(num_iters) :: errors
    real(ccs_real), dimension(num_iters) :: errors_upwind
    real(ccs_real), dimension(num_iters) :: errors_luds
    real(ccs_real), dimension(num_iters) :: errors_cd
    real(ccs_real), dimension(num_iters) :: errors_gamma

    real(ccs_real), dimension(num_iters) :: refinements
    real(ccs_real) :: order

    real(ccs_real), dimension(3, 2) :: rvecs, grads

    type(advection_kernel) :: advection
    type(upwind_advection_kernel) :: upwind_advection
    type(luds_advection_kernel) :: luds_advection
    type(cd_advection_kernel) :: cd_advection
    type(gamma_advection_kernel) :: gamma_advection


    ! Initialising values
    x_P = [0, 0, 0]
    phi_P = phi(x_P(1))
    interpol_factor = 2.0_ccs_real / 3.0_ccs_real
    rvecs = 0.0_ccs_real
    grads = 0.0_ccs_real

    ! Refinement loop
    do i = 1, num_iters
        dx = 1.0_ccs_real / real(i, ccs_real)**2
        x_N = x_P + [dx, 0.0_ccs_real, 0.0_ccs_real]
        x_f = x_P + [dx/3, 0.0_ccs_real, 0.0_ccs_real]

        phi_N = phi(dx)
        u_f = interpol_factor*u(x_P(1)) + (1-interpol_factor)*u(x_N(1))

        coeffs = advection %eval_coeffs(u_f)
        rhs = advection%eval_explicit(u_f, interpol_factor, rvecs, grads)

        ! Compute discretisation error
        call get_error(coeffs, rhs, x_P, x_N, x_f, errors(i))

        ! Upwind Kernel
        coeffs = upwind_advection %eval_coeffs(u_f)
        rhs = upwind_advection%eval_explicit(u_f, interpol_factor, rvecs, grads)
        call get_error(coeffs, rhs, x_P, x_N, x_f, errors_upwind(i))

        ! Linearise Upwind Kernel
        coeffs = luds_advection %eval_coeffs(u_f)
        rhs = luds_advection%eval_explicit(u_f, interpol_factor, rvecs, grads)
        call get_error(coeffs, rhs, x_P, x_N, x_f, errors_luds(i))

        ! Central Difference Kernel
        coeffs = cd_advection %eval_coeffs(u_f)
        rhs = cd_advection%eval_explicit(u_f, interpol_factor, rvecs, grads)
        call get_error(coeffs, rhs, x_P, x_N, x_f, errors_cd(i))

        ! Gamma Kernel
        coeffs = gamma_advection %eval_coeffs(u_f)
        rhs = gamma_advection%eval_explicit(u_f, interpol_factor, rvecs, grads)
        call get_error(coeffs, rhs, x_P, x_N, x_f, errors_gamma(i))


        refinements(i) = dx
    end do

    ! Compute convergence order
    call compute_order(refinements, errors, order)
    call assert_gt(order, advection%get_order() * 0.95_ccs_real, "Convergence order not preserved by advection kernel")

    call compute_order(refinements, errors_upwind, order)
    call assert_gt(order, upwind_advection%get_order() * 0.95_ccs_real, "Convergence order not preserved by upwind advection kernel")

    call compute_order(refinements, errors_luds, order)
    call assert_gt(order, luds_advection%get_order() * 0.95_ccs_real, "Convergence order not preserved by linearised upwind advection kernel")

    call compute_order(refinements, errors_cd, order)
    call assert_gt(order, cd_advection%get_order() * 0.95_ccs_real, "Convergence order not preserved by central difference advection kernel")

    call compute_order(refinements, errors_gamma, order)
    call assert_gt(order, gamma_advection%get_order() * 0.95_ccs_real, "Convergence order not preserved by gamma advection kernel")

contains

    function u(x)
        ! Velocity field
        real(ccs_real) :: u
        real(ccs_real), intent(in) :: x

        u = cos(x)
    end function u

    function phi(x)
        ! Function being tested
        real(ccs_real) :: phi
        real(ccs_real), intent(in) :: x

        phi = sin(2*x + 17)
    end function phi

    subroutine get_error(coeffs, rhs, x_P, x_N, x_f, error)
        ! Computes discretisation error by comparing with analytical value 
        real(ccs_real), intent(in) :: coeffs(2), rhs
        real(ccs_real), intent(in) :: x_P(3)
        real(ccs_real), intent(in) :: x_N(3)
        real(ccs_real), intent(in) :: x_f(3)
        real(ccs_real), intent(out) :: error
        real(ccs_real) :: analytical

        analytical = u(x_f(1)) * phi(x_f(1))
        error = abs(analytical - (coeffs(1) * phi(x_P(1)) + coeffs(2) * phi(x_N(1)) + rhs))
    end subroutine get_error

end program test_advection_kernel
