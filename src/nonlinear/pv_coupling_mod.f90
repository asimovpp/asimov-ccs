!> @brief Module file pv_coupling.mod
!
!> @details An interface to pressure-velocity coupling methods (SIMPLE, etc)

module pv_coupling

    use kinds, only : accs_real, accs_int

    implicit none

    private

    public :: solve_nonlinear

    interface

    module subroutine solve_nonlinear(u, v, p)
        type(field), intent(inout) :: u, v, p
    end subroutine solve_nonlinear

    end interface

end module pv_coupling