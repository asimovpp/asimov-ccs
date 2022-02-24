!> @brief Submodule file pv_coupling_simple.smod
!
!> @details Implementation of the SIMPLE algorithm for pressure-velocity coupling.

submodule (pv_coupling) pv_coupling_simple

  use fv, only : compute_fluxes

  implicit none

  contains

  module subroutine solve_nonlinear(u, v, p, M, solution, source, cell_mesh)

    ! Arguments
    class(field), intent(inout) :: u, v, p
    class(matrix), intent(inout) :: M
    class(vector), intent(inout) :: solution, source
    type(mesh), intent(in) :: cell_mesh

    ! Local variables
    integer(accs_int) :: i

    outerloop: do i = it_start, it_end
      call calculate_momentum(u, v, cell_mesh, bcs, cps, M, source)

      call calculate_pressure_correction()

      call update_pressure()

      ! Todo:
      !call calculate_scalars()

      call check_convergence()
      if (converged) then
        exit outerloop
      endif

    end do outerloop

  end subroutine

  subroutine calculate_momentum(u, v, cell_mesh, bcs, cps, M, vec)

    ! Arguments
    type(field), intent(inout)    :: u, v
    type(mesh), intent(in)        :: cell_mesh
    type(bc_config), intent(in)   :: bcs
    integer(accs_int), intent(in) :: cps
    class(matrix), intent(inout)  :: M
    class(vector), intent(inout)  :: vec

    ! Local variables
    
    ! Calculate fluxes
    call compute_fluxes(u, v, cell_mesh, bcs, cps, M, vec)

    ! Assembly of coefficient matrix and source vector
    call update(M)
    call update(vec)

    ! Create linear solver
    call set_linear_system(scalar_linear_system, vec, scalar, M, par_env)
    call create_solver(scalar_linear_system, scalar_solver)

    ! Solve the linear system
    call solve(scalar_solver)

  end subroutine calculate_momentum


  subroutine calculate_pressure_correction()

    


  end subroutine calculate_pressure_correction


  subroutine update_pressure()


  end subroutine update_pressure


  subroutine calculate_scalars()


  end subroutine calculate_scalars


  subroutine check_convergence(converged)

    ! Arguments
    logical, intent(inout) :: converged

  end subroutine check_convergence

end submodule pv_coupling_simple