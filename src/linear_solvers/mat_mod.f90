!> @brief Module file mat_mod.f90
!
!> @details Provides the interface to matrix objects.
module mat

  use kinds, only : ccs_int, ccs_real
  use types, only : ccs_matrix, matrix_spec, matrix_values, ccs_mesh, ccs_vector
  use parallel_types, only: parallel_environment

  implicit none

  private

  public :: create_matrix
  public :: finalise_matrix
  public :: update_matrix
  public :: begin_update_matrix
  public :: end_update_matrix
  public :: set_matrix_values
  public :: set_eqn
  public :: pack_one_matrix_coefficient
  public :: initialise_matrix
  public :: set_matrix_size
  public :: set_nnz
  public :: mat_axpy
  public :: mat_norm
  public :: get_matrix_diagonal
  public :: set_matrix_diagonal
  public :: zero_matrix
  
  interface

     !> @brief Interface to create a new matrix object.
     !
     !> @param[in]  mat_dat - contains information about how the matrix should be allocated
     !> @param[out] M       - the matrix object
     module subroutine create_matrix(mat_dat, M)
       type(matrix_spec), intent(in) :: mat_dat
       class(ccs_matrix), allocatable, intent(out) :: M
     end subroutine

    module subroutine finalise_matrix(M)
      class(ccs_matrix), intent(inout) :: M
    end subroutine

     !> @brief Interface to perform a parallel update of a matrix.
     !
     !> @param[in/out] M - the matrix
     module subroutine update_matrix(M)
       class(ccs_matrix), intent(inout) :: M
     end subroutine

     !> @brief Interface to begin a parallel update of a matrix.
     !
     !> @param[in/out] M - the matrix
     !
     !> @details Begins the parallel update to allow overlapping comms and compute.
     module subroutine begin_update_matrix(M)
       class(ccs_matrix), intent(inout) :: M
     end subroutine

     !> @brief Interface to end a parallel update of a matrix.
     !
     !> @param[in/out] M - the matrix
     !
     !> @details Ends the parallel update to allow overlapping comms and compute.
     module subroutine end_update_matrix(M)
       class(ccs_matrix), intent(inout) :: M
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
     module subroutine pack_one_matrix_coefficient(row_entry, col_entry, row, col, coeff, mat_coeffs)
       integer(ccs_int), intent(in) :: row_entry
       integer(ccs_int), intent(in) :: col_entry
       integer(ccs_int), intent(in) :: row
       integer(ccs_int), intent(in) :: col
       real(ccs_real), intent(in) :: coeff
       type(matrix_values), intent(inout) :: mat_coeffs
     end subroutine pack_one_matrix_coefficient

    !> @brief Interface to perform the AXPY matrix operation.
    !
    !> @details Performs the AXPY operation
    !!          y[i] = a * x[i] + y[i]
    !
    !> @param[in]     alpha - a scalar value
    !> @param[in]     x     - an input matrix
    !> @param[in,out] y     - matrix serving as input, overwritten with result
    module subroutine mat_axpy(alpha, x, y)
      real(ccs_real), intent(in) :: alpha
      class(ccs_matrix), intent(in) :: x
      class(ccs_matrix), intent(inout) :: y
    end subroutine

    !> @brief Interface to compute the norm of a matrix
    !
    !> @param[in]  m         - the matrix
    !> @param[in]  norm_type - which norm to compute? Currently supported is the 2 norm:
    !!                         norm_type=2.
    !> @param[out] n         - the computed norm returned as the result of the function
    !!                         call.
    module function mat_norm(M, norm_type) result(n)
      class(ccs_matrix), intent(in) :: M
      integer(ccs_int), intent(in) :: norm_type
      real(ccs_real) :: n
    end function
     
     !> @brief Interface to set values in a matrix.
     !
     !> @param[in]     mat_values - contains the values, their indices and the mode to use when setting
     !!                             them.
     !> @param[in/out] M          - the matrix
     module subroutine set_matrix_values(mat_values, M)
       type(matrix_values), intent(in) :: mat_values
       class(ccs_matrix), intent(inout) :: M
     end subroutine

     !> @brief Interface to set equation
     !
     !> @param[in]  rows - array of (global) row indices to set the equation on
     !> @param[in/out] M - the matrix
     !
     !> @details Sets equations in a system of equations by zeroing out the corresponding row in the
     !!          system matrix and setting the diagonal to one such that the solution is given by
     !!          the corresponding entry in the right-hand side vector.
     module subroutine set_eqn(rows, M)
       integer(ccs_int), dimension(:), intent(in) :: rows
       class(ccs_matrix), intent(inout) :: M
     end subroutine

    !> @brief Constructor for default matrix values
    !
    !> param[in/out] matrix_descriptor - the initialised matrix values
     module subroutine initialise_matrix(matrix_descriptor)
      type(matrix_spec), intent(inout) :: matrix_descriptor
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
    end subroutine

    !> @brief Setter for matrix number of non-zeros
    !
    !> param[in] nnz                   - the number of non-zeros
    !> param[in/out] matrix_descriptor - the matrix data object
    module subroutine set_nnz(nnz, matrix_descriptor)
      integer(ccs_int), intent(in) :: nnz
      type(matrix_spec), intent(inout) :: matrix_descriptor
    end subroutine

    !> @brief Extract matrix diagonal elements into a vector
    !
    !> @param[in]  M - the matrix
    !> @param[out] D - a vector containing the diagonal elements of M
    module subroutine get_matrix_diagonal(M, D)
      class(ccs_matrix), intent(in)  :: M
      class(ccs_vector), intent(inout) :: D
    end subroutine get_matrix_diagonal

    module subroutine set_matrix_diagonal(D, M)
      class(ccs_vector), intent(in) :: D
      class(ccs_matrix), intent(inout) :: M
    end subroutine set_matrix_diagonal

    module subroutine zero_matrix(M)
      class(ccs_matrix), intent(inout) :: M
    end subroutine zero_matrix
    
  end interface
  
end module mat
