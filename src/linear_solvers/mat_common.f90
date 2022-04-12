submodule (mat) mat_common

  implicit none

contains

  !> @brief Constructor for default matrix values
  !
  !> param[in/out] mat_properties - the initialised matrix values
  module subroutine initialise_matrix(mat_properties)
    type(matrix_spec), intent(inout) :: mat_properties
    mat_properties%par_env => null()
    mat_properties%par_env => null()
  end subroutine initialise_matrix

  !> @brief Setter for global matrix size
  !
  !> param[in] par_env                - the parallel environment where 
  !!                                    the matrix resides
  !> param[in] geometry               - the mesh object
  !> param[in/out] mat_properties  - the matrix data object
  module subroutine set_matrix_size(par_env, geometry, mat_properties)
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    class(ccs_mesh), target, intent(in) :: geometry
    type(matrix_spec), intent(inout) :: mat_properties

    mat_properties%mesh => geometry
    mat_properties%par_env => par_env
  end subroutine set_matrix_size

  !> @brief Setter for matrix number of non-zeros
  !
  !> param[in] nnz                   - the number of non-zeros
  !> param[in/out] mat_properties - the matrix data object
  module subroutine set_nnz(nnz, mat_properties)
    integer(ccs_int), intent(in) :: nnz
    type(matrix_spec), intent(inout) :: mat_properties

    mat_properties%nnz = nnz
  end subroutine  set_nnz

end submodule
