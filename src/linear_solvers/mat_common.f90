submodule (mat) mat_common

  implicit none

contains

  !> @brief Constructor for default matrix values
  !
  !> param[in/out] matrix_descriptor - the initialised matrix values
  module subroutine initialise_matrix(matrix_descriptor)
    type(matrix_init_data), intent(inout) :: matrix_descriptor
    matrix_descriptor%rglob = -1
    matrix_descriptor%cglob = -1
    matrix_descriptor%rloc = -1
    matrix_descriptor%cloc = -1
    matrix_descriptor%nnz = -1
    matrix_descriptor%par_env => null()
  end subroutine initialise_matrix

  !> @brief Setter for global matrix size
  !
  !> param[in/out] matrix_descriptor  - the matrix data object
  !> param[in] rows                   - the global number of rows
  !> param[in] columns                - the global number of columns
  !> param[in] par_env                - the parallel environment where 
  !!                                    the matrix resides
  module subroutine set_global_matrix_size(matrix_descriptor, rows, columns, par_env)
    type(matrix_init_data), intent(inout) :: matrix_descriptor
    integer(accs_int), intent(in) :: rows
    integer(accs_int), intent(in) :: columns
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    matrix_descriptor%rglob = rows
    matrix_descriptor%cglob = columns
    matrix_descriptor%rloc = -1
    matrix_descriptor%cloc = -1
    matrix_descriptor%par_env => par_env
  end subroutine

  !> @brief Setter for local matrix size
  !
  !> param[in/out] matrix_descriptor  - the matrix data object
  !> param[in] rows                   - the local number of rows
  !> param[in] columns                - the local number of columns
  !> param[in] par_env                - the parallel environment where
  !!                                    the matrix resides
  module subroutine set_local_matrix_size(matrix_descriptor, rows, columns, par_env)
    type(matrix_init_data), intent(inout) :: matrix_descriptor
    integer(accs_int), intent(in) :: rows
    integer(accs_int), intent(in) :: columns
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    matrix_descriptor%rloc = rows
    matrix_descriptor%cloc = columns
    matrix_descriptor%rglob = -1
    matrix_descriptor%cglob = -1
    matrix_descriptor%par_env => par_env
  end subroutine

  !> @brief Setter for matrix number of non-zeros
  !
  !> param[in/out] matrix_descriptor - the matrix data object
  !> param[in] nnz                   - the number of non-zeros
  !> param[in] par_env               - the parallel environment where
  !!                                   the matrix resides
  module subroutine set_nnz(matrix_descriptor, nnz)
    type(matrix_init_data), intent(inout) :: matrix_descriptor
    integer(accs_int), intent(in) :: nnz
    matrix_descriptor%nnz = nnz
  end subroutine  

  module procedure create_matrix_values
    allocate(val_dat%rglob(nrows))
    allocate(val_dat%cglob(nrows))
  end procedure create_matrix_values

  module procedure set_matrix_values_mode
    val_dat%mode = mode
  end procedure set_matrix_values_mode
  
  module subroutine set_matrix_values_entry(val, val_dat)

    use constants, only : add_mode, insert_mode

    real(accs_real), intent(in) :: val
    type(matrix_values), intent(inout) :: val_dat
    
    associate(x => val_dat%val(val_dat%current_entry), &
         mode => val_dat%mode)
      if (mode == insert_mode) then
        x = val
      else if (mode == add_mode) then
        x = x + val
      else
        print *, "ERROR: Unrecognised entry mode ", mode
        stop
      end if
    end associate
    
  end subroutine set_matrix_values_entry
  
end submodule 