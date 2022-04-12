submodule (mat) mat_common

  implicit none

contains

  !> @brief Constructor for default matrix values
  !
  !> param[in/out] matrix_descriptor - the initialised matrix values
  module subroutine initialise_matrix(matrix_descriptor)
    type(matrix_spec), intent(inout) :: matrix_descriptor
    matrix_descriptor%par_env => null()
    matrix_descriptor%par_env => null()
  end subroutine initialise_matrix

  !> @brief Setter for global matrix size
  !
  !> param[in] par_env                - the parallel environment where 
  !!                                    the matrix resides
  !> param[in] mesh               - the mesh object
  !> param[in/out] matrix_descriptor  - the matrix data object
  module subroutine set_matrix_size(par_env, mesh, matrix_descriptor)
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    class(ccs_mesh), target, intent(in) :: mesh
    type(matrix_spec), intent(inout) :: matrix_descriptor

    matrix_descriptor%mesh => mesh
    matrix_descriptor%par_env => par_env
  end subroutine set_matrix_size

  !> @brief Setter for matrix number of non-zeros
  !
  !> param[in] nnz                   - the number of non-zeros
  !> param[in/out] matrix_descriptor - the matrix data object
  module subroutine set_nnz(nnz, matrix_descriptor)
    integer(ccs_int), intent(in) :: nnz
    type(matrix_spec), intent(inout) :: matrix_descriptor

    matrix_descriptor%nnz = nnz
  end subroutine  set_nnz

end submodule
