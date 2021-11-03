submodule (mat) mat_common

  use types, only: matrix_init_data

  implicit none

contains

  !> @brief Constructor for default matrix values
  !
  !> param[in/out] mat     - the initialised matrix values
  module subroutine initialise_matrix(mat)
    type(matrix_init_data), intent(inout) :: mat
    mat%rglob = -1
    mat%cglob = -1
    mat%rloc = -1
    mat%cloc = -1
    mat%nnz = -1
    mat%par_env => null()
  end subroutine initialise_matrix

end submodule