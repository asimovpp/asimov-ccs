!>  Module file solver.mod
!
!>  An interface to linear solver objects.

module solver

  use types, only : linear_solver, equation_system, ccs_vector, ccs_matrix
  use parallel_types, only: parallel_environment
  use vec, only : vec_axpy, vec_norm
  use mat, only : mat_axpy, mat_norm
  
  implicit none

  private

  public :: create_solver
  public :: solve
  public :: initialise_equation_system
  public :: set_equation_system
  public :: axpy
  public :: norm
  public :: set_solver_method
  public :: set_solver_precon
  
  interface

    !>  Interface to create a new solver object.
    !
    !> @param[in]  equation_system linear_system   - Data structure containing equation system to be solved.
    !> @param[out] linear_solver solver - The linear solver returned allocated.
    module subroutine create_solver(linear_system, solver)
      type(equation_system), intent(in) :: linear_system
      class(linear_solver), allocatable, intent(out) :: solver
    end subroutine

    !>  Interface to solve the linear system in a solver.
    !
    !> @param[in/out] linear_solver solver - The linear solver object.
    module subroutine solve(solver)
      class(linear_solver), intent(inout) :: solver
    end subroutine

    !>  Constructor for default linear system
    module subroutine initialise_equation_system(lin_sys)
      type(equation_system), intent(inout) :: lin_sys
    end subroutine initialise_equation_system

    !>  Setter for the linear system
    !
    !> @param[in] par_env       - the parallel environment where the linear 
    !!                            system resides
    !> @param[in] rhs           - the right hand side vector
    !> @param[in] solution      - the solution vector
    !> @param[in] mat           - the matrix
    !> @param[in/out] lin_sys   - the linear system
    module subroutine set_equation_system(par_env, rhs, solution, mat, lin_sys)
      class(parallel_environment), allocatable, target, intent(in) :: par_env
      class(ccs_vector), allocatable, target, intent(in) :: rhs
      class(ccs_vector), allocatable, target, intent(in) :: solution
      class(ccs_matrix), allocatable, target, intent(in) :: mat
      type(equation_system), intent(inout) :: lin_sys
    end subroutine

    !> Interface to set the primary method of a linear solver
    module subroutine set_solver_method(method_name, solver)
      character(len=*), intent(in) :: method_name       !> String naming the linear solver to be used.
      class(linear_solver), intent(inout) :: solver !> The linear solver object
    end subroutine set_solver_method

    !> Interface to set the preconditioner of a linear solver
    module subroutine set_solver_precon(precon_name, solver)
      character(len=*), intent(in) :: precon_name       !> String naming the preconditioner to be used.
      class(linear_solver), intent(inout) :: solver !> The linear solver object
    end subroutine set_solver_precon
    
  end interface
  
  !>  Generic interface to perform the AXPY operation (a*x + y)
  interface axpy
    module procedure vec_axpy
    module procedure mat_axpy
  end interface axpy
  
  !>  Generic interface to compute the norm of an element
  interface norm
    module procedure vec_norm
    module procedure mat_norm
  end interface norm
  
end module solver
