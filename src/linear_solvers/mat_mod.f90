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
  public :: update_matrix
  public :: begin_update_matrix
  public :: end_update_matrix
  public :: set_matrix_values
  public :: set_eqn
  public :: pack_one_matrix_coefficient

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

     module subroutine pack_one_matrix_coefficient(mat_coeffs, row_entry, col_entry, row, col, coeff)
       type(matrix_values), intent(inout) :: mat_coeffs
       integer(accs_int), intent(in) :: row_entry
       integer(accs_int), intent(in) :: col_entry
       integer(accs_int), intent(in) :: row
       integer(accs_int), intent(in) :: col
       real(accs_real), intent(in) :: coeff
     end subroutine pack_one_matrix_coefficient
     
     !> @brief Interface to set values in a matrix.
     !
     !> @param[in]     mat_values - contains the values, their indices and the mode to use when setting
     !!                             them.
     !> @param[in/out] M          - the matrix
     module subroutine set_matrix_values(mat_values, M)
       type(matrix_values), intent(in) :: mat_values
       class(matrix), intent(inout) :: M
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
       integer(accs_int), dimension(:), intent(in) :: rows
       class(matrix), intent(inout) :: M
     end subroutine
     
  end interface
  
end module mat
