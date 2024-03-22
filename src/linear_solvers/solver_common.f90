submodule(solver) solver_common

  implicit none

contains

  !> Constructor for default linear system
  pure module subroutine initialise_equation_system(lin_sys)
    type(equation_system), intent(inout) :: lin_sys

    lin_sys%solution => null()
    lin_sys%rhs => null()
    lin_sys%matrix => null()
    lin_sys%par_env => null()
  end subroutine initialise_equation_system

  !> Setter for the linear system
  module subroutine set_equation_system(par_env, rhs, solution, mat, lin_sys, name)
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< the parallel environment where the linear system resides
    class(ccs_vector), allocatable, target, intent(in) :: rhs               !< the right hand side vector
    class(ccs_vector), allocatable, target, intent(in) :: solution          !< the solution vector
    class(ccs_matrix), allocatable, target, intent(in) :: mat               !< the matrix
    type(equation_system), intent(inout) :: lin_sys                         !< the linear system
    character(len=*), optional, intent(in) :: name                          !< name of the equation system

    lin_sys%rhs => rhs
    lin_sys%solution => solution
    lin_sys%matrix => mat
    lin_sys%par_env => par_env
    if (present(name)) lin_sys%name = name
  end subroutine

end submodule
