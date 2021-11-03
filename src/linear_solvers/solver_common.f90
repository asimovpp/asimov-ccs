submodule (solver) solver_common

  use types, only: linear_system

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

  end submodule