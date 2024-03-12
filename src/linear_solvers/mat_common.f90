submodule(mat) mat_common
#include "ccs_macros.inc"

  use utils, only: exit_print, str
  use error_codes

  implicit none

contains

  !> Set number of rows in working set specifier
  pure module subroutine set_matrix_values_spec_nrows(nrows, val_spec)

    ! Arguments
    integer(ccs_int), intent(in) :: nrows               !< Number of rows to be used in the working set
    type(matrix_values_spec), intent(inout) :: val_spec !< The current working set specifier object

    val_spec%nrows = nrows
  end subroutine set_matrix_values_spec_nrows

  !> Set number of columns in working set specifier
  pure module subroutine set_matrix_values_spec_ncols(ncols, val_spec)

    ! Arguments
    integer(ccs_int), intent(in) :: ncols               !< Number of columns to be used in the working set
    type(matrix_values_spec), intent(inout) :: val_spec !< The current working setspecifier object

    val_spec%ncols = ncols
  end subroutine set_matrix_values_spec_ncols

  !> Constructor for default matrix values
  pure module subroutine initialise_matrix(mat_properties)

    ! Arguments
    type(matrix_spec), intent(inout) :: mat_properties !< The initialised matrix values

    mat_properties%par_env => null()
    mat_properties%par_env => null()
  end subroutine initialise_matrix

  !> Setter for global matrix size
  module subroutine set_matrix_size(par_env, mesh, mat_properties)

    ! Arguments
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< the parallel environment where the matrix resides
    class(ccs_mesh), target, intent(in) :: mesh                             !< the mesh object
    type(matrix_spec), intent(inout) :: mat_properties                      !< the matrix data object

    mat_properties%mesh => mesh
    mat_properties%par_env => par_env
  end subroutine set_matrix_size

  !> Setter for matrix number of non-zeros
  pure module subroutine set_nnz(nnz, mat_properties)

    ! Arguments
    integer(ccs_int), intent(in) :: nnz                !< the number of non-zeros
    type(matrix_spec), intent(inout) :: mat_properties !< the matrix data object

    mat_properties%nnz = nnz
  end subroutine set_nnz

  !> Constructor for matrix values object
  pure module subroutine create_matrix_values(val_spec, val_dat)

    ! Arguments
    type(matrix_values_spec), intent(in) :: val_spec !< Object describing the size (nrows, ncol) of working set.
    type(matrix_values), intent(out) :: val_dat      !< The working set object.

    associate (nrows => val_spec%nrows, ncols => val_spec%ncols)
      allocate (val_dat%global_row_indices(nrows))
      allocate (val_dat%global_col_indices(ncols))
      allocate (val_dat%values(nrows * ncols))
    end associate

    ! Run clear entries to ensure correct default values.
    call clear_matrix_values_entries(val_dat)
  end subroutine create_matrix_values

  !> Set the storage mode.
  pure module subroutine set_matrix_values_mode(mode, val_dat)
    integer(ccs_int), intent(in) :: mode          !< The storage mode
    type(matrix_values), intent(inout) :: val_dat !< The object storing the working set

    val_dat%setter_mode = mode
  end subroutine set_matrix_values_mode

  !v Store a coefficient in the current working set at the current row,col coordinate, using the
  !  current storage mode.
  pure module subroutine set_matrix_values_entry(val, val_dat)

    use constants, only: add_mode, insert_mode

    ! Arguments
    real(ccs_real), intent(in) :: val             !< The coefficient value
    type(matrix_values), intent(inout) :: val_dat !< The object storing the working set

    ! Local
    integer(ccs_int) :: current_entry ! Logically 2D index of current row,col coordinate.
    ! XXX: This may be PETSc-specific, but seems sensible for now

    ! Locate row,col coordinate in logically-2D array
    associate (row => val_dat%current_row, &
               col => val_dat%current_col, &
               nrows => size(val_dat%global_row_indices), &
               ncols => size(val_dat%global_col_indices))
      current_entry = (row - 1) * ncols + col
    end associate

    ! Store value according to working set mode.
    associate (x => val_dat%values(current_entry), &
               mode => val_dat%setter_mode)
      if (mode == insert_mode) then
        x = val
      else if (mode == add_mode) then
        x = x + val
      else
        error stop unknown_mode ! Unrecognised entry mode
      end if

    end associate

  end subroutine set_matrix_values_entry

end submodule
