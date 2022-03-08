!> @brief Module file pv_coupling.mod
!
!> @details An interface to pressure-velocity coupling methods (SIMPLE, etc)

module pv_coupling

    use kinds, only : accs_real, accs_int
    use types, only: vector, matrix, field, mesh

    implicit none

    private

    public :: solve_nonlinear

    interface

    module subroutine solve_nonlinear(u, v, p, pp, M, solution, source, cell_mesh)
        class(field), intent(inout) :: u, v, p, pp
        class(matrix), intent(inout) :: M
        class(vector), intent(inout) :: solution, source
        type(mesh), intent(in) :: cell_mesh
    end subroutine solve_nonlinear

    end interface

end module pv_coupling