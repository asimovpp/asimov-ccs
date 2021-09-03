!> @brief Module file accs_solver.mod
!>
!> @details An interface to linear solver objects.

module accs_solver

  use accs_types, only : linear_solver, linear_system
  
  implicit none

  private

  public :: create_solver, solve
  
  interface

     !> @brief Interface to create a new solver object.
     !>
     !> @param[in] linear_system eqsys   - Data structure containing equation system to be solved.
     !> @param[out] linear_solver solver - The linear solver returned allocated.
     module subroutine create_solver(eqsys, solver)
       type(linear_system), intent(in) :: eqsys
       class(linear_solver), allocatable, intent(out) :: solver
     end subroutine

     !> @brief Interface to solve the linear system in a solver.
     !>
     !> @param[in/out] linear_solver solver - The linear solver object.
     module subroutine solve(solver)
       class(linear_solver), intent(inout) :: solver
     end subroutine

  end interface
  
end module accs_solver
