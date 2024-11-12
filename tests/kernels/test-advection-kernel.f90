
program test_advection_kernel
    !! Test the advection kernel by running a refinement loop, computing the discretisation
    !!  error (difference between analytical and computed values), and checking the order of
    !!  convergence is the theoretical order for the scheme
    implicit none

    use kinds, only: ccs_real


    integer(ccs_int) :: i, num_iters=10
    real(ccs_real) :: phi_P, interpol_factor
    real(ccs_real) :: phi_N
    real(ccs_real) :: dx
    real(ccs_real), dimension(3) :: x_N, x_P, x_f
    real(ccs_real), dimension(2) :: coeffs
    real(ccs_real) :: rhs
    real(ccs_real) :: u_f
    real(ccs_real), dimension(num_iters) :: errors
    real(ccs_real), dimension(num_iters) :: refinements
    real(ccs_real) :: order


    type(advection_kernel) :: advection

    ! Initialising values
    x_P = [0, 0, 0]
    phi_P = phi(x_P(1))
    interpol_factor = 2_ccs_real / 3_ccs_real

    ! Refinement loop
    do i = 1, num_iters
        dx = 1.0_ccs_real / real(i)**2
        x_N = x_P + [dx, 0, 0]
        x_f = x_P + [dx/3, 0, 0]

        phi_N = phi(dx)
        u_f = interpol_factor*u(x_P(1)) + (1-interpol_factor)*u(x_N(1))

        call advection%eval([phi_P, phi_N], u_f, x_P, x_N, x_f, &
                            is_boundary=.false., coeffs, rhs)

        ! Compute discretisation error
        call get_error(coeffs, rhs, x_P, x_N, x_f, errors(i))
        refinements(i) = dx
    end do

    ! Compute convergence order
    call get_order(refinements, errors, order)

    call assert_gt(order, advection%get_order()*0.95, "Convergence order not preserved by advection kernel")

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
        real(ccs_real) :: analytical, error

        analytical = u(x_f(1)) * phi(x_f(1))
        error = abs(analytical - (coeffs(1) * phi(x_P(1)) + coeffs(2) * phi(x_N(1)) + rhs))
    end subroutine get_error

end program test_advection_kernel