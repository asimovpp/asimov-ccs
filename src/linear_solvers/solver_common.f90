submodule (solver) solver_common

  implicit none

contains

  !> @brief Constructor for default linear system
  module subroutine initialise_linear_system(lin_sys)
    type(linear_system), intent(inout) :: lin_sys
    
    lin_sys%sol => null()
    lin_sys%rhs => null()
    lin_sys%M => null()
    lin_sys%par_env => null()
  end subroutine initialise_linear_system

  !> @brief Setter for the linear system
  !
  !> @param[in/out] lin_sys   - the linear system
  !> @param[in] rhs           - the right hand side vector
  !> @param[in] solution      - the solution vector
  !> @param[in] mat           - the matrix
  !> @param[in] par_env       - the parallel environment where the linear 
  !!                            system resides
  module subroutine set_linear_system(lin_sys, rhs, solution, mat, par_env)
    type(linear_system), intent(inout) :: lin_sys
    class(vector), allocatable, target, intent(in) :: rhs
    class(vector), allocatable, target, intent(in) :: solution
    class(matrix), allocatable, target, intent(in) :: mat
    class(parallel_environment), allocatable, target, intent(in) :: par_env

    lin_sys%rhs => rhs
    lin_sys%sol => solution
    lin_sys%M => mat
    lin_sys%par_env => par_env
  end subroutine

end submodule