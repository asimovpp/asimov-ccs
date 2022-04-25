submodule (mat) mat_common

  implicit none

contains

  !>  Constructor for default matrix values
  module subroutine initialise_matrix(mat_properties)
    type(matrix_spec), intent(inout) :: mat_properties  !< the initialised matrix values
    mat_properties%par_env => null()
    mat_properties%par_env => null()
  end subroutine initialise_matrix

  !>  Setter for global matrix size
  module subroutine set_matrix_size(par_env, mesh, mat_properties)
    class(parallel_environment), allocatable, target, intent(in) :: par_env   !< the parallel environment where the matrix resides
    class(ccs_mesh), target, intent(in) :: mesh                               !< the mesh object
    type(matrix_spec), intent(inout) :: mat_properties                        !< the matrix data object

    mat_properties%mesh => mesh
    mat_properties%par_env => par_env
  end subroutine set_matrix_size

  !>  Setter for matrix number of non-zeros
  module subroutine set_nnz(nnz, mat_properties)
    integer(ccs_int), intent(in) :: nnz                 !< the number of non-zeros
    type(matrix_spec), intent(inout) :: mat_properties  !< the matrix data object

    mat_properties%nnz = nnz
  end subroutine  set_nnz

end submodule
