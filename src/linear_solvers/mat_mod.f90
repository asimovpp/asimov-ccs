!v Module file mat_mod.f90
!
!  Provides the interface to matrix objects.
module mat

  use kinds, only: ccs_int, ccs_real
  use types, only: ccs_matrix, matrix_spec, matrix_values, matrix_values_spec, ccs_mesh, ccs_vector
  use parallel_types, only: parallel_environment

  implicit none

  private

  public :: create_matrix
  public :: finalise_matrix
  public :: get_info_matrix
  public :: set_matrix_values
  public :: clear_matrix_values_entries
  public :: set_matrix_values_entry
  public :: create_matrix_values
  public :: set_matrix_values_mode
  public :: set_matrix_values_row
  public :: set_matrix_values_col
  public :: set_matrix_values_spec_nrows
  public :: set_matrix_values_spec_ncols
  public :: update_matrix
  public :: begin_update_matrix
  public :: end_update_matrix
  public :: set_eqn
  public :: initialise_matrix
  public :: set_matrix_size
  public :: set_nnz
  public :: mat_axpy
  public :: mat_norm
  public :: mat_vec_product
  public :: get_matrix_diagonal
  public :: set_matrix_diagonal
  public :: add_matrix_diagonal
  public :: zero_matrix

  interface

    !> Interface to create a new matrix object.
    module subroutine create_matrix(mat_properties, M, name)
      type(matrix_spec), intent(in) :: mat_properties  !< contains information about
      !< how the matrix should be allocated
      class(ccs_matrix), allocatable, intent(out) :: M !< the matrix object
      character(len=*), optional, intent(in) :: name !< name of the matrix object
    end subroutine

    module subroutine finalise_matrix(M)
      class(ccs_matrix), intent(inout) :: M
    end subroutine

    !> Returns information about matrix storage (number of nonzeros, memory, etc.)
    module subroutine get_info_matrix(M)
      class(ccs_matrix), intent(inout) :: M
    end subroutine

    !> Interface to set values in a matrix.
    module subroutine set_matrix_values(mat_values, M)
      type(matrix_values), intent(in) :: mat_values !< contains the values, their indices
      !< and the mode to use when setting them.
      class(ccs_matrix), intent(inout) :: M         !< the matrix
    end subroutine

    !> Clear working set of values to begin new working set.
    pure module subroutine clear_matrix_values_entries(val_dat)

      ! Arguments
      type(matrix_values), intent(inout) :: val_dat !< Working set object

    end subroutine clear_matrix_values_entries

    !v Store a coefficient in the current working set at the current row,col coordinate, using the
    !  current storage mode.
    pure module subroutine set_matrix_values_entry(val, val_dat)

      ! Arguments
      real(ccs_real), intent(in) :: val             !< The coefficient value
      type(matrix_values), intent(inout) :: val_dat !< The object storing the working set

    end subroutine set_matrix_values_entry

    !> Set the storage mode.
    pure module subroutine set_matrix_values_mode(mode, val_dat)

      ! Arguments
      integer(ccs_int), intent(in) :: mode          !< The storage mode
      type(matrix_values), intent(inout) :: val_dat !< The object storing the working set

    end subroutine set_matrix_values_mode

    !v Interface to set the row currently being worked on by matrix values.
    !
    !  Sets the current row in the matrix value object, the implementation of this is
    !  backend-dependent as it should immediately convert to the correct indexing
    !  (whether that's 0, 1 or X-based) as used by the backend.
    pure module subroutine set_matrix_values_row(row, val_dat)
      integer(ccs_int), intent(in) :: row           !< the row
      type(matrix_values), intent(inout) :: val_dat !< the matrix values object
    end subroutine set_matrix_values_row

    pure module subroutine set_matrix_values_col(col, val_dat)
      integer(ccs_int), intent(in) :: col
      type(matrix_values), intent(inout) :: val_dat
    end subroutine set_matrix_values_col

    !> Set number of rows in working set specifier
    pure module subroutine set_matrix_values_spec_nrows(nrows, val_spec)

      ! Arguments
      integer(ccs_int), intent(in) :: nrows               !< Number of rows to be used in the working set
      type(matrix_values_spec), intent(inout) :: val_spec !< The current working set specifier object

    end subroutine set_matrix_values_spec_nrows

    !> Set number of columns in working set specifier
    pure module subroutine set_matrix_values_spec_ncols(ncols, val_spec)

      ! Arguments
      integer(ccs_int), intent(in) :: ncols               !< Number of columns to be used in the working set
      type(matrix_values_spec), intent(inout) :: val_spec !< The current working setspecifier object

    end subroutine set_matrix_values_spec_ncols

    !> Interface to create a matrix values object.
    pure module subroutine create_matrix_values(val_spec, val_dat)
      type(matrix_values_spec), intent(in) :: val_spec !< how many rows will be set?
      type(matrix_values), intent(out) :: val_dat      !< the matrix values object
    end subroutine create_matrix_values

    !> Interface to perform a parallel update of a matrix.
    module subroutine update_matrix(M)
      class(ccs_matrix), intent(inout) :: M !< the matrix
    end subroutine

    !v Interface to begin a parallel update of a matrix.
    !
    !  Begins the parallel update to allow overlapping comms and compute.
    module subroutine begin_update_matrix(M)
      class(ccs_matrix), intent(inout) :: M !< the matrix
    end subroutine

    !v Interface to end a parallel update of a matrix.
    !
    !  Ends the parallel update to allow overlapping comms and compute.
    module subroutine end_update_matrix(M)
      class(ccs_matrix), intent(inout) :: M !< the matrix
    end subroutine

    !v Interface to perform the AXPY matrix operation.
    !
    !  Performs the AXPY operation
    !
    !       y[i] = a * x[i] + y[i]
    module subroutine mat_axpy(alpha, x, y)
      real(ccs_real), intent(in) :: alpha   !< a scalar value
      class(ccs_matrix), intent(in) :: x    !< an input matrix
      class(ccs_matrix), intent(inout) :: y !< matrix serving as input, overwritten with result
    end subroutine

    !> Interface to compute the norm of a matrix
    module function mat_norm(M, norm_type) result(n)
      class(ccs_matrix), intent(in) :: M        !< the matrix
      integer(ccs_int), intent(in) :: norm_type !< which norm to compute?
      !< Currently supported is the 2 norm: norm_type=2.
      real(ccs_real) :: n !< the computed norm returned as the result of the function call.
    end function

    !> Compute matrix-vector product
    module subroutine mat_vec_product(M, x, y)
      class(ccs_matrix), intent(in) :: M    !< the input matrix
      class(ccs_vector), intent(in) :: x    !< the input vector
      class(ccs_vector), intent(inout) :: y !< the output (result) vector
    end subroutine

    !> Interface to set equation
    !
    !  Sets equations in a system of equations by zeroing out the corresponding row in the
    !  system matrix and setting the diagonal to one such that the solution is given by
    !  the corresponding entry in the right-hand side vector.
    module subroutine set_eqn(global_rows, M)
      integer(ccs_int), dimension(:), intent(in) :: global_rows !< array of (global) row indices to set the equation on
      class(ccs_matrix), intent(inout) :: M !< the matrix
    end subroutine

    !> Constructor for default matrix values
    pure module subroutine initialise_matrix(mat_properties)
      type(matrix_spec), intent(inout) :: mat_properties !< the initialised matrix values
    end subroutine initialise_matrix

    !> Setter for global matrix size
    module subroutine set_matrix_size(par_env, mesh, mat_properties)
      class(parallel_environment), allocatable, target, intent(in) :: par_env !< the parallel environment where the matrix resides
      class(ccs_mesh), target, intent(in) :: mesh !< the mesh object
      type(matrix_spec), intent(inout) :: mat_properties !< the matrix data object
    end subroutine

    !> Setter for matrix number of non-zeros
    pure module subroutine set_nnz(nnz, mat_properties)
      integer(ccs_int), intent(in) :: nnz                !< the number of non-zeros
      type(matrix_spec), intent(inout) :: mat_properties !< the matrix data object
    end subroutine

    !> Extract matrix diagonal elements into a vector
    module subroutine get_matrix_diagonal(M, D)
      class(ccs_matrix), intent(in) :: M    !< the matrix
      class(ccs_vector), intent(inout) :: D !< a vector containing the diagonal elements of M
    end subroutine get_matrix_diagonal

    module subroutine set_matrix_diagonal(D, M)
      class(ccs_vector), intent(in) :: D
      class(ccs_matrix), intent(inout) :: M
    end subroutine set_matrix_diagonal

    !> Add a vector to the matrix diagonal
    module subroutine add_matrix_diagonal(D, M)
      class(ccs_vector), intent(in) :: D      !< the vector containing matrix diagonal elements
      class(ccs_matrix), intent(inout) :: M   !< the matrix
    end subroutine add_matrix_diagonal
    
    module subroutine zero_matrix(M)
      class(ccs_matrix), intent(inout) :: M
    end subroutine zero_matrix
  end interface

end module mat
