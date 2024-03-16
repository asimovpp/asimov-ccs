!v Module file solver.mod
!
!  An interface to linear solver objects.

module solver

  use types, only: linear_solver, equation_system, ccs_vector, ccs_matrix
  use parallel_types, only: parallel_environment
  use vec, only: vec_axpy, vec_norm
  use mat, only: mat_axpy, mat_norm

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

    !> Interface to create a new solver object.
    module subroutine create_solver(linear_system, solver)
      type(equation_system), intent(in) :: linear_system !< Data structure containing
      !< equation system to be solved.
      class(linear_solver), allocatable, intent(inout) :: solver !< The linear solver returned allocated.
    end subroutine

    !> Interface to solve the linear system in a solver.
    module subroutine solve(solver)
      class(linear_solver), intent(inout) :: solver !< The linear solver object.
    end subroutine

    !> Constructor for default linear system
    pure module subroutine initialise_equation_system(lin_sys)
      type(equation_system), intent(inout) :: lin_sys
    end subroutine initialise_equation_system

    !> Setter for the linear system
    module subroutine set_equation_system(par_env, rhs, solution, mat, lin_sys, name)
      class(parallel_environment), allocatable, target, intent(in) :: par_env !< the parallel environment
      !< where the linear system resides
      class(ccs_vector), allocatable, target, intent(in) :: rhs               !< the right hand side vector
      class(ccs_vector), allocatable, target, intent(in) :: solution          !< the solution vector
      class(ccs_matrix), allocatable, target, intent(in) :: mat               !< the matrix
      type(equation_system), intent(inout) :: lin_sys                         !< the linear system
      character(len=*), optional, intent(in) :: name                          !< name of the equation system
    end subroutine

    !> Interface to set the primary method of a linear solver
    module subroutine set_solver_method(method_name, solver)
      character(len=*), intent(in) :: method_name   !< String naming the linear solver to be used.
      class(linear_solver), intent(inout) :: solver !< The linear solver object
    end subroutine set_solver_method

    !> Interface to set the preconditioner of a linear solver
    module subroutine set_solver_precon(precon_name, solver)
      character(len=*), intent(in) :: precon_name   !< String naming the preconditioner to be used.
      class(linear_solver), intent(inout) :: solver !< The linear solver object
    end subroutine set_solver_precon

  end interface

  !> Generic interface to perform the AXPY operation (a*x + y)
  interface axpy
    module procedure vec_axpy
    module procedure mat_axpy
  end interface axpy

  !> Generic interface to compute the norm of an element
  interface norm
    module procedure vec_norm
    module procedure mat_norm
  end interface norm

end module solver
