!> @brief Module file mat_mod.f90
!
!> @details Provides the interface to matrix objects.
module mat

  use kinds, only : accs_int, accs_real
  use types, only : matrix, matrix_init_data, matrix_values
  use parallel_types, only: parallel_environment

  implicit none

  private

  public :: create_matrix
  public :: finalise_matrix
  public :: set_matrix_values
  public :: clear_matrix_values_entries
  public :: set_matrix_values_entry
  public :: create_matrix_values
  public :: set_matrix_values_mode
  public :: set_matrix_values_row
  public :: update_matrix
  public :: begin_update_matrix
  public :: end_update_matrix
  public :: set_eqn
  public :: pack_one_matrix_coefficient
  public :: initialise_matrix
  public :: set_global_matrix_size
  public :: set_local_matrix_size
  public :: set_nnz

  interface

     !> @brief Interface to create a new matrix object.
     !
     !> @param[in]  mat_dat - contains information about how the matrix should be allocated
     !> @param[out] M       - the matrix object
     module subroutine create_matrix(mat_dat, M)
       type(matrix_init_data), intent(in) :: mat_dat
       class(matrix), allocatable, intent(out) :: M
     end subroutine

    module subroutine finalise_matrix(M)
      class(matrix), intent(inout) :: M
    end subroutine

     !> @brief Interface to set values in a matrix.
     !
     !> @param[in]     mat_values - contains the values, their indices and the mode to use when setting
     !!                             them.
     !> @param[in/out] M          - the matrix
    module subroutine set_matrix_values(mat_values, M)
      type(matrix_values), intent(in) :: mat_values
      class(matrix), intent(inout) :: M
    end subroutine
    
    module subroutine clear_matrix_values_entries(val_dat)
      type(matrix_values), intent(inout) ::val_dat
    end subroutine clear_matrix_values_entries

    module subroutine set_matrix_values_entry(val, val_dat)
      real(accs_real), intent(in) :: val
      type(matrix_values), intent(inout) :: val_dat
    end subroutine set_matrix_values_entry

    module subroutine set_matrix_values_mode(mode, val_dat)
      integer(accs_int), intent(in) :: mode
      type(matrix_values), intent(inout) :: val_dat
    end subroutine set_matrix_values_mode

    !> @brief Interface to set the row currently being worked on by matrix values.
    !
    !> @description Sets the current row in the maitrx value object, the implementation of this is
    !!              backend-dependent as it should immediately convert to the correct indexing
    !!              (whether that's 0, 1 or X-based) as used by the backend.
    !
    !> @param[in]     row     - the row
    !> @param[in,out] val_dat - the matrix values object
    module subroutine set_matrix_values_row(row, val_dat)
      integer(accs_int), intent(in) :: row
      type(matrix_values), intent(inout) :: val_dat
    end subroutine set_matrix_values_row

    !> @brief Interface to create a maitrx values object.
    !
    !> @param[in]  nrows   - how many rows will be set?
    !> @param[out] val_dat - the matrix values object
    module subroutine create_matrix_values(nrows, val_dat)
      integer(accs_int), intent(in) :: nrows
      type(matrix_values), intent(out) :: val_dat
    end subroutine create_matrix_values

    !> @brief Interface to perform a parallel update of a matrix.
     !
     !> @param[in/out] M - the matrix
     module subroutine update_matrix(M)
       class(matrix), intent(inout) :: M
     end subroutine

     !> @brief Interface to begin a parallel update of a matrix.
     !
     !> @param[in/out] M - the matrix
     !
     !> @details Begins the parallel update to allow overlapping comms and compute.
     module subroutine begin_update_matrix(M)
       class(matrix), intent(inout) :: M
     end subroutine

     !> @brief Interface to end a parallel update of a matrix.
     !
     !> @param[in/out] M - the matrix
     !
     !> @details Ends the parallel update to allow overlapping comms and compute.
     module subroutine end_update_matrix(M)
       class(matrix), intent(inout) :: M
     end subroutine

     !> @brief Interface to store one matrix coefficient and its index for later setting.
     !
     !> @param[in/out] mat_coeffs - object storing the coefficients, their indices and mode to use
     !!                             when setting them.
     !> @param[in]     row_entry  - which entry in the row indices to set?
     !> @param[in]     col_entry  - which entry in the column indices to set?
     !> @param[in]     row        - matrix row index
     !> @param[in]     col        - matrix column index
     !> @param[in]     coeff      - matrix coefficient
     !
     !> @details Stores a matrix coefficient and associated row and column indices for later
     !!          setting, ensuring they are set appropriately for the backend.
     module subroutine pack_one_matrix_coefficient(mat_coeffs, row_entry, col_entry, row, col, coeff)
      type(matrix_values), intent(inout) :: mat_coeffs
      integer(accs_int), intent(in) :: row_entry
      integer(accs_int), intent(in) :: col_entry
      integer(accs_int), intent(in) :: row
      integer(accs_int), intent(in) :: col
      real(accs_real), intent(in) :: coeff
    end subroutine pack_one_matrix_coefficient

    !> @brief Interface to set equation
     !
     !> @param[in]  rows - array of (global) row indices to set the equation on
     !> @param[in/out] M - the matrix
     !
     !> @details Sets equations in a system of equations by zeroing out the corresponding row in the
     !!          system matrix and setting the diagonal to one such that the solution is given by
     !!          the corresponding entry in the right-hand side vector.
     module subroutine set_eqn(rows, M)
       integer(accs_int), dimension(:), intent(in) :: rows
       class(matrix), intent(inout) :: M
     end subroutine

    !> @brief Constructor for default matrix values
    !
    !> param[in/out] matrix_descriptor - the initialised matrix values
     module subroutine initialise_matrix(matrix_descriptor)
      type(matrix_init_data), intent(inout) :: matrix_descriptor
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
    end subroutine

  end interface
  
end module mat
