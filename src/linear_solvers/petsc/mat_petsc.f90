submodule(mat) mat_petsc
#include "ccs_macros.inc"

  use kinds, only: ccs_err
  use petsctypes, only: matrix_petsc, vector_petsc
  use parallel_types_mpi, only: parallel_environment_mpi
  use petscmat, only: MatAssemblyBegin, MatAssemblyEnd, MAT_FLUSH_ASSEMBLY
  use petsc, only: ADD_VALUES, INSERT_VALUES
  use utils, only: debug_print, str, update, exit_print
  use error_codes

  implicit none

contains

  !> Create a new PETSc matrix object.
  module subroutine create_matrix(mat_properties, M, name)

    use mpi

    use petsc, only: PETSC_DETERMINE, PETSC_NULL_INTEGER
    use petscmat, only: MatCreate, MatSetSizes, MatSetFromOptions, MatSetUp, &
                        MatSeqAIJSetPreallocation, MatMPIAIJSetPreallocation

    use meshing, only: get_local_num_cells

    type(matrix_spec), intent(in) :: mat_properties   !< contains information about how the matrix should be allocated
    class(ccs_matrix), allocatable, intent(out) :: M  !< the matrix object
    character(len=*), optional, intent(in) :: name    !< name of the matrix object

    integer(ccs_int) :: local_num_cells

    integer(ccs_err) :: ierr  ! Error code

    allocate (matrix_petsc :: M)

    select type (M)
    type is (matrix_petsc)

      M%modeset = .false.
      if (present(name)) M%name = name

      select type (par_env => mat_properties%par_env)
      type is (parallel_environment_mpi)

        call MatCreate(par_env%comm, M%M, ierr)

        associate (mesh => mat_properties%mesh)
          call get_local_num_cells(local_num_cells)
          call MatSetSizes(M%M, local_num_cells, local_num_cells, &
                           PETSC_DETERMINE, PETSC_DETERMINE, ierr)
        end associate

        if (ierr == 0) then
          M%allocated = .true.
        end if

        if (present(name)) then
          call MatSetOptionsPrefix(M%M, M%name // ':', ierr)
        end if
        call MatSetFromOptions(M%M, ierr)

        if (mat_properties%nnz < 1) then
          if (par_env%proc_id == par_env%root) then
            call dprint("WARNING: No matrix preallocation set, potentially inefficient.")
          end if
          call MatSetUp(M%M, ierr)
        else
          call MatSeqAIJSetPreallocation(M%M, mat_properties%nnz, PETSC_NULL_INTEGER, ierr)
          call MatMPIAIJSetPreallocation(M%M, mat_properties%nnz, PETSC_NULL_INTEGER, mat_properties%nnz - 1, &
                                         PETSC_NULL_INTEGER, ierr)
        end if

      class default
        call error_abort("Unknown parallel environment")

      end select

    class default
      call error_abort("Unsupported matrix type")

    end select

  end subroutine

  module subroutine finalise_matrix(M)

    use petscmat, only: MAT_FINAL_ASSEMBLY

    class(ccs_matrix), intent(inout) :: M

    integer(ccs_err) :: ierr

    select type (M)
    type is (matrix_petsc)
      call MatAssemblyBegin(M%M, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(M%M, MAT_FINAL_ASSEMBLY, ierr)
    end select

  end subroutine finalise_matrix

  !> Returns information about matrix storage (number of nonzeros, memory, etc.) 
  ! see https://petsc.org/release/manualpages/Mat/MatInfo/ for all the available fields
  module subroutine get_info_matrix(M)

    use petscmat, only: MAT_INFO_SIZE, MatGetInfo, MAT_INFO_MEMORY, MAT_INFO_NZ_ALLOCATED, MAT_LOCAL, &
       MAT_INFO_NZ_USED, MAT_INFO_NZ_UNNEEDED

    class(ccs_matrix), intent(inout) :: M
    double precision, dimension(MAT_INFO_SIZE) :: info

    integer(ccs_err) :: ierr

    select type (M)
    type is (matrix_petsc)
      call MatGetInfo(M%M, MAT_LOCAL, info, ierr)
      print *, "---"
      print *, "nnz allocated: ", info(MAT_INFO_NZ_ALLOCATED)
      print *, "nnz used: ", info(MAT_INFO_NZ_USED)
      print *, "nnz unneeded: ", info(MAT_INFO_NZ_UNNEEDED)
    end select

  end subroutine get_info_matrix

  !> Perform a parallel update of a PETSc matrix.
  module subroutine update_matrix(M)

    class(ccs_matrix), intent(inout) :: M   !< the matrix

    select type (M)
    type is (matrix_petsc)

      call begin_update_matrix(M)
      call end_update_matrix(M)

    class default
      call error_abort("Unsupported matrix type")

    end select

  end subroutine

  !v Begin a parallel update of a PETSc matrix.
  !
  !  Begins the parallel update to allow overlapping comms and compute.
  module subroutine begin_update_matrix(M)

    class(ccs_matrix), intent(inout) :: M   !< the matrix

    integer(ccs_err) :: ierr ! Error code

    select type (M)
    type is (matrix_petsc)

      call MatAssemblyBegin(M%M, MAT_FLUSH_ASSEMBLY, ierr)

    class default
      call error_abort("Unsupported matrix type")

    end select

  end subroutine

  !v End a parallel update of a PETSc matrix.
  !
  !  Ends the parallel update to allow overlapping comms and compute.
  module subroutine end_update_matrix(M)

    class(ccs_matrix), intent(inout) :: M   !< the matrix

    integer(ccs_err) :: ierr ! Error code

    select type (M)
    type is (matrix_petsc)

      call MatAssemblyEnd(M%M, MAT_FLUSH_ASSEMBLY, ierr)

      M%modeset = .false. ! It's safe to change modes now
    class default
      call error_abort("Unsupported matrix type")

    end select

  end subroutine

  !> Set values in a PETSc matrix.
  module subroutine set_matrix_values(mat_values, M)

    use petscmat, only: MatSetValues
    use constants, only: insert_mode, add_mode

    type(matrix_values), intent(in) :: mat_values   !< contains the values, their indices and the mode to use when setting them.
    class(ccs_matrix), intent(inout) :: M           !< the matrix

    integer(ccs_int) :: nrows, ncols ! number of rows/columns
    integer(ccs_int) :: mode ! Add or insert values?

    integer(ccs_err) :: ierr ! Error code

    associate (ridx => mat_values%global_row_indices, &
               cidx => mat_values%global_col_indices, &
               val => mat_values%values, &
               matmode => mat_values%setter_mode)

      select type (M)
      type is (matrix_petsc)

        if (M%modeset) then
          if (matmode /= M%mode) then
            call error_abort("ERROR: changing matrix mode without updating")
          end if
        else
          M%mode = matmode
          M%modeset = .true.
        end if

        nrows = size(ridx)
        ncols = size(cidx)
        if (nrows * ncols /= size(val)) then
          call error_abort("Invalid matrix values.")
        end if
        if (matmode == add_mode) then
          mode = ADD_VALUES
        else if (matmode == insert_mode) then
          mode = INSERT_VALUES
        else
          call error_abort("Unknown mode.")
        end if

        call MatSetValues(M%M, nrows, ridx, ncols, cidx, val, mode, ierr)

      class default
        call error_abort("Unknown matrix type.")

      end select

    end associate

  end subroutine

  !v Set equation
  !
  !  Sets equations in a system of equations by zeroing out the corresponding row in the
  !  system matrix and setting the diagonal to one such that the solution is given by
  !  the corresponding entry in the right-hand side vector.  module subroutine set_eqn(rows, M)
  module subroutine set_eqn(global_rows, M)

    use petsc, only: PETSC_NULL_VEC
    use petscmat, only: MatZeroRows

    integer(ccs_int), dimension(:), intent(in) :: global_rows  !< array of (global) row indices to set the equation on
    class(ccs_matrix), intent(inout) :: M                      !< the matrix

    integer(ccs_err) :: ierr

    select type (M)
    type is (matrix_petsc)

      call MatZeroRows(M%M, size(global_rows), global_rows, 1.0_ccs_real, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)

    class default
      call error_abort("Unknown matrix type.")

    end select

  end subroutine

  !> Clear working set of values to begin new working set.
  pure module subroutine clear_matrix_values_entries(val_dat)
    type(matrix_values), intent(inout) :: val_dat !< Working set object

    val_dat%global_row_indices(:) = -1 ! PETSc ignores -ve indices, used as "empty" indicator
    val_dat%global_col_indices(:) = -1 ! PETSc ignores -ve indices, used as "empty" indicator
    val_dat%values(:) = 0.0_ccs_real
  end subroutine clear_matrix_values_entries

  !> Set working row.
  pure module subroutine set_matrix_values_row(row, val_dat)

    ! Arguments
    integer(ccs_int), intent(in) :: row           !< Which (global) row to work on?
    type(matrix_values), intent(inout) :: val_dat !< Object recording values and their coordinates

    ! Local
    integer(ccs_int), dimension(1) :: rglobs ! Temporary array mapping rows to indices in the
    ! current working set. N.B. the dimension of this
    ! array must match the rank of
    ! matrix_values%global_row_indices!
    integer(ccs_int) :: i         ! The mapped index in the current working set
    logical :: new_entry          ! Flag to indicate if we are revisiting a row
    integer(ccs_int) :: petsc_row ! The (zero-indexed) row as used by PETSc

    petsc_row = row - 1 ! PETSc is zero-indexed
    new_entry = .false.

    rglobs = findloc(val_dat%global_row_indices, petsc_row, kind=ccs_int)
    i = rglobs(1) ! We want the first entry
    if (i == 0) then
      new_entry = .true.
    end if

    if (new_entry) then
      rglobs = findloc(val_dat%global_row_indices, -1_ccs_int, kind=ccs_int)
      i = rglobs(1) ! We want the first entry
      if (i == 0) then
        error stop no_free_entry ! Couldn't find a free entry in matrix values
      end if
    end if

    val_dat%current_row = i
    val_dat%global_row_indices(i) = petsc_row

  end subroutine set_matrix_values_row

  !> Set working column.
  pure module subroutine set_matrix_values_col(col, val_dat)

    ! Arguments
    integer(ccs_int), intent(in) :: col           !< Which (global) column to work on ?
    type(matrix_values), intent(inout) :: val_dat !< Object recording values and their coordinates

    ! Local
    integer(ccs_int), dimension(1) :: cglobs ! Temporary array mapping columns to indices in the
    ! current working set. N.B. the dimension of this
    ! array must match the rank of
    ! matrix_values%global_col_indices!
    integer(ccs_int) :: i         ! The mapped index in the current working set
    logical :: new_entry          ! Flag to indicate if we are revisiting a column
    integer(ccs_int) :: petsc_col ! The (zero-indexed) column as used by PETSc

    petsc_col = col - 1 ! PETSc is zero-indexed
    new_entry = .false.

    cglobs = findloc(val_dat%global_col_indices, petsc_col, kind=ccs_int)
    i = cglobs(1) ! We want the first entry
    if (i == 0) then
      new_entry = .true.
    end if

    if (new_entry) then
      cglobs = findloc(val_dat%global_col_indices, -1_ccs_int, kind=ccs_int)
      i = cglobs(1) ! We want the first entry
      if (i == 0) then
        error stop no_free_entry ! Couldn't find a free column entry in matrix values
      end if
    end if

    val_dat%current_col = i
    val_dat%global_col_indices(i) = petsc_col

  end subroutine set_matrix_values_col

  !v Perform the AXPY matrix operation using PETSc
  !
  !         y[i] = alpha * x[i] + y[i]
  module subroutine mat_axpy(alpha, x, y)

    use petscmat, only: MatAXPY, DIFFERENT_NONZERO_PATTERN

    real(ccs_real), intent(in) :: alpha     !< a scalar value
    class(ccs_matrix), intent(in) :: x      !< a PETSc input matrix
    class(ccs_matrix), intent(inout) :: y   !< PETSc matrix serving as input, overwritten with result

    integer(ccs_err) :: ierr ! Error code

    select type (x)
    type is (matrix_petsc)

      select type (y)
      type is (matrix_petsc)

        ! PETSc performs AXPY as YPAX, with result stored in Y.
        call MatAXPY(y%M, alpha, x%M, DIFFERENT_NONZERO_PATTERN, ierr)

      class default
        call error_abort("Unknown matrix type.")

      end select

    class default
      call error_abort("Unknown matrix type.")

    end select

  end subroutine

  !> Compute the norm of a PETSc matrix
  module function mat_norm(M, norm_type) result(n)

    use petscmat, only: NORM_1, NORM_FROBENIUS, NORM_INFINITY, MatNorm

    class(ccs_matrix), intent(in) :: M         !< the PETSc matrix
    integer(ccs_int), intent(in) :: norm_type  !< which norm to compute

    real(ccs_real) :: n      !< The computed norm
    integer(ccs_err) :: ierr ! Error code

    n = 0.0_ccs_real ! initialise norm to 0

    select type (M)
    type is (matrix_petsc)

      if (norm_type == 1) then
        call MatNorm(M%M, NORM_1, n, ierr)
      else if (norm_type == 2) then
        call MatNorm(M%M, NORM_FROBENIUS, n, ierr)
      else if (norm_type == 3) then
        call MatNorm(M%M, NORM_INFINITY, n, ierr)
      else
        call error_abort("ERROR: unknown matrix norm type " // str(norm_type))
      end if

    class default
      call error_abort("Type unhandled")
    end select

  end function

  !> Extract the diagonal elements of a matrix and store in a vector
  module subroutine get_matrix_diagonal(M, D)

    use petscmat, only: MatGetDiagonal

    class(ccs_matrix), intent(in) :: M     !< the PETSc matrix
    class(ccs_vector), intent(inout) :: D  !< the PETSc vector containing matrix diagonal elements

    integer(ccs_err) :: ierr ! Error code

    select type (M)
    type is (matrix_petsc)

      select type (D)
      type is (vector_petsc)
        call MatGetDiagonal(M%M, D%v, ierr)

      class default
        call error_abort("Unknown vector type.")
      end select

    class default
      call error_abort("Unknown matrix type.")
    end select

  end subroutine

  !> Store a vector in the matrix diagonal
  module subroutine set_matrix_diagonal(D, M)
    use petscmat, only: MatDiagonalSet

    class(ccs_vector), intent(in) :: D      !< the PETSc vector containing matrix diagonal elements
    class(ccs_matrix), intent(inout) :: M   !< the PETSc matrix

    integer(ccs_err) :: ierr

    select type (M)
    type is (matrix_petsc)

      select type (D)
      type is (vector_petsc)
        call MatDiagonalSet(M%M, D%v, INSERT_VALUES, ierr)

      class default
        call error_abort("Unknown vector type.")
      end select

    class default
      call error_abort("Unknown matrix type.")
    end select

  end subroutine set_matrix_diagonal

  !> Add a vector to the matrix diagonal
  module subroutine add_matrix_diagonal(D, M)
    use petscmat, only: MatDiagonalSet

    class(ccs_vector), intent(in) :: D      !< the PETSc vector containing matrix diagonal elements
    class(ccs_matrix), intent(inout) :: M   !< the PETSc matrix

    integer(ccs_err) :: ierr

    select type (M)
    type is (matrix_petsc)

      select type (D)
      type is (vector_petsc)
        call MatDiagonalSet(M%M, D%v, ADD_VALUES, ierr)

      class default
        call error_abort("Unknown vector type.")
      end select

    class default
      call error_abort("Unknown matrix type.")
    end select

  end subroutine add_matrix_diagonal

  !> Overwite a matrix with zeros.
  module subroutine zero_matrix(M)

    use petscmat, only: MatZeroEntries

    class(ccs_matrix), intent(inout) :: M   !< the PETSc matrix

    integer(ccs_err) :: ierr ! Error code

    select type (M)
    type is (matrix_petsc)
      call MatZeroEntries(M%M, ierr)
    class default
      call error_abort("Unknown matrix type.")

    end select

  end subroutine zero_matrix

  !> Compute matrix-vector product
  module subroutine mat_vec_product(M, x, y)

    use petscmat, only: MatMult

    ! Arguments
    class(ccs_matrix), intent(in) :: M
    class(ccs_vector), intent(in) :: x
    class(ccs_vector), intent(inout) :: y

    ! Local variables
    integer(ccs_err) :: ierr

    select type (M)
    type is (matrix_petsc)

      select type (x)
      type is (vector_petsc)

        select type (y)
        type is (vector_petsc)

          call MatMult(M%M, x%v, y%v, ierr)

        class default
          call error_abort("Unknown vector type.")
        end select

      class default
        call error_abort("Unknown vector type.")
      end select

    class default
      call error_abort("Unknown matrix type.")
    end select

  end subroutine mat_vec_product

end submodule mat_petsc
