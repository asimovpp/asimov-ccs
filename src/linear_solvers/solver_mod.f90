!> @brief Module file solver.mod
!
!> @details An interface to linear solver objects.

module solver

  use types, only : linear_solver, linear_system
  use parallel_types, only: parallel_environment
  
  implicit none

  private

  public :: create_solver
  public :: solve
  
  interface

     !> @brief Interface to create a new solver object.
     !
     !> @param[in]  linear_system eqsys   - Data structure containing equation system to be solved.
     !> @param[out] linear_solver solver - The linear solver returned allocated.
     module subroutine create_solver(eqsys, solver)
       type(linear_system), intent(in) :: eqsys
       class(linear_solver), allocatable, intent(out) :: solver
     end subroutine

     !> @brief Interface to solve the linear system in a solver.
     !
     !> @param[in/out] linear_solver solver - The linear solver object.
     module subroutine solve(solver)
       class(linear_solver), intent(inout) :: solver
     end subroutine

  end interface
  
end module solver
