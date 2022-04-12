submodule (solver) solver_common

  implicit none

contains

  !> @brief Constructor for default linear system
  module subroutine initialise_equation_system(lin_sys)
    type(equation_system), intent(inout) :: lin_sys
    
    lin_sys%solution => null()
    lin_sys%rhs => null()
    lin_sys%matrix => null()
    lin_sys%par_env => null()
  end subroutine initialise_equation_system

  !> @brief Setter for the linear system
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

    lin_sys%rhs => rhs
    lin_sys%solution => solution
    lin_sys%matrix => mat
    lin_sys%par_env => par_env
  end subroutine

end submodule
