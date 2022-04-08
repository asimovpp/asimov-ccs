!> @brief Module file solver.mod
!
!> @details An interface to linear solver objects.

module solver

  use types, only : linear_solver, linear_system, vector, matrix
  use parallel_types, only: parallel_environment
  use vec, only : vec_axpy, vec_norm
  use mat, only : mat_axpy, mat_norm
  
  implicit none

  private

  public :: create_solver
  public :: solve
  public :: initialise_linear_system
  public :: set_linear_system
  public :: axpy
  public :: norm
  
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

    !> @brief Constructor for default linear system
    module subroutine initialise_linear_system(lin_sys)
      type(linear_system), intent(inout) :: lin_sys
    end subroutine initialise_linear_system

    !> @brief Setter for the linear system
    !
    !> @param[in] par_env       - the parallel environment where the linear 
    !!                            system resides
    !> @param[in] rhs           - the right hand side vector
    !> @param[in] solution      - the solution vector
    !> @param[in] mat           - the matrix
    !> @param[in/out] lin_sys   - the linear system
    module subroutine set_linear_system(par_env, rhs, solution, mat, lin_sys)
      class(parallel_environment), allocatable, target, intent(in) :: par_env
      class(vector), allocatable, target, intent(in) :: rhs
      class(vector), allocatable, target, intent(in) :: solution
      class(matrix), allocatable, target, intent(in) :: mat
      type(linear_system), intent(inout) :: lin_sys
    end subroutine

  end interface
  
  !> @brief Generic interface to perform the AXPY operation (a*x + y)
  interface axpy
    module procedure vec_axpy
    module procedure mat_axpy
  end interface axpy
  
  !> @brief Generic interface to compute the norm of an element
  interface norm
    module procedure vec_norm
    module procedure mat_norm
  end interface norm
  
end module solver
