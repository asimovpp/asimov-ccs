!> @brief Module file pv_coupling.mod
!
!> @details An interface to pressure-velocity coupling methods (SIMPLE, etc)

module pv_coupling

    use kinds, only : accs_int
    use types, only: field, mesh
    use parallel_types, only: parallel_environment

    implicit none

    private

    public :: solve_nonlinear

    interface

    module subroutine solve_nonlinear(par_env, cell_mesh, it_start, it_end, u, v, p, pp)
        class(parallel_environment), intent(in) :: par_env
        type(mesh), intent(in) :: cell_mesh
        integer(accs_int), intent(in) :: it_start, it_end
        class(field), intent(inout) :: u, v, p, pp
    end subroutine solve_nonlinear

    end interface

end module pv_coupling