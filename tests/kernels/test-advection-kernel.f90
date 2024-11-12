
program test_advection_kernel
    implicit none

    use kinds, only: ccs_real


    integer(ccs_int) :: i, num_iters=10
    real(ccs_real) :: phi_P, phi_P_prime, interpol_factor
    real(ccs_real) :: phi_N
    real(ccs_real) :: dx
    real(ccs_real), dimension(3) :: N_loc, P_loc, f_loc
    real(ccs_real), dimension(2) :: coeffs
    real(ccs_real) :: rhs
    real(ccs_real), dimension(num_iters) :: errors
    real(ccs_real), dimension(num_iters) :: refinements
    real(ccs_real) :: order


    type(advection_kernel) :: advection


    phi_P = phi(0)
    phi_P_prime = grad_phi(0)
    interpol_factor = 2_ccs_real / 3_ccs_real

    do i = 1, num_iters
        dx = 1.0_ccs_real / real(i)**2
        N_loc = P_loc + [dx, 0, 0]
        f_loc = P_loc + [dx/3, 0, 0]

        phi_N = phi(dx)

        call advection%eval([phi_P, phi_N], [u_P, u_N], P_loc, N_loc, f_loc, &
                            interpol_factor, is_boundary=.false., coeffs, rhs)

        call get_error(coeffs, rhs, f_loc, errors(i))
        refinements(i) = dx
    end do

    call get_order(refinements, errors, order)

    call assert_gt(order, advection%get_order()*0.95, "Convergence order not preserved by advection kernel")

contains

    function u(x)
        real(ccs_real) :: u
        real(ccs_real), intent(in) :: x

        u = cos(x)
    end function u

    function phi(x)
        real(ccs_real) :: phi
        real(ccs_real), intent(in) :: x

        phi = sin(2*x + 17)
    end function phi

    function grad_phi(x)
        real(ccs_real) :: grad_phi
        real(ccs_real), intent(in) :: x

        grad_phi = 2*cos(2*x + 17)
    end function grad_phi

    subroutine get_error(coeffs, rhs, f_loc, error)
        real(ccs_real), intent(in) :: coeffs(2), rhs
        real(ccs_real), intent(in) :: f_loc(3)
        real(ccs_real), intent(out) :: error
        real(ccs_real) :: analytical, error

        analytical = u(f_loc(1)) * grad_phi(f_loc(1))
        error = abs(analytical - (rhs - coeffs(2) * phi(f_loc(1))) / coeffs(1))
    end subroutine get_error

end program test_advection_kernel