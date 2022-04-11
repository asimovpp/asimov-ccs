submodule (mat) mat_common

  implicit none

contains

  !> @brief Constructor for default matrix values
  !
  !> param[in/out] matrix_descriptor - the initialised matrix values
  module subroutine initialise_matrix(matrix_descriptor)
    type(matrix_init_data), intent(inout) :: matrix_descriptor
    matrix_descriptor%par_env => null()
    matrix_descriptor%par_env => null()
  end subroutine initialise_matrix

  !> @brief Setter for global matrix size
  !
  !> param[in] par_env                - the parallel environment where 
  !!                                    the matrix resides
  !> param[in] geometry               - the mesh object
  !> param[in/out] matrix_descriptor  - the matrix data object
  module subroutine set_matrix_size(par_env, geometry, matrix_descriptor)
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    class(mesh), target, intent(in) :: geometry
    type(matrix_init_data), intent(inout) :: matrix_descriptor

    matrix_descriptor%mesh => geometry
    matrix_descriptor%par_env => par_env
  end subroutine set_matrix_size

  !> @brief Setter for matrix number of non-zeros
  !
  !> param[in] nnz                   - the number of non-zeros
  !> param[in/out] matrix_descriptor - the matrix data object
  module subroutine set_nnz(nnz, matrix_descriptor)
    integer(ccs_int), intent(in) :: nnz
    type(matrix_init_data), intent(inout) :: matrix_descriptor

    matrix_descriptor%nnz = nnz
  end subroutine  set_nnz

end submodule
