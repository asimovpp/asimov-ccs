submodule (mat) mat_common

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

  module subroutine set_global_matrix_size(mat, rows, columns, par_env)
    type(matrix_init_data), intent(inout) :: mat
    integer(accs_int), intent(in) :: rows
    integer(accs_int), intent(in) :: columns
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    mat%rglob = rows
    mat%cglob = columns
    mat%par_env => par_env
  end subroutine

  module subroutine set_local_matrix_size(mat, rows, columns, par_env)
    type(matrix_init_data), intent(inout) :: mat
    integer(accs_int), intent(in) :: rows
    integer(accs_int), intent(in) :: columns
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    mat%rloc = rows
    mat%cloc = columns
    mat%par_env => par_env
  end subroutine

  module subroutine set_nnz(mat, nnz)
    type(matrix_init_data), intent(inout) :: mat
    integer(accs_int), intent(in) :: nnz
    mat%nnz = nnz
  end subroutine  

end submodule