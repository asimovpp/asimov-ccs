submodule (mat) mat_petsc

  use kinds, only : ccs_err
  use petsctypes, only : matrix_petsc, vector_petsc
  use parallel_types_mpi, only: parallel_environment_mpi
  use petscmat, only: MatAssemblyBegin, MatAssemblyEnd, MAT_FLUSH_ASSEMBLY
  use petsc, only : ADD_VALUES, INSERT_VALUES
  
  implicit none

contains

  !> @brief Create a new PETSc matrix object.
  !
  !> @param[in]  mat_dat - contains information about how the matrix should be allocated
  !> @param[out] M       - the matrix object
  module subroutine create_matrix(mat_dat, M)

    use mpi
    
    use petsc, only : PETSC_DETERMINE, PETSC_NULL_INTEGER
    use petscmat, only : MatCreate, MatSetSizes, MatSetFromOptions, MatSetUp, &
                         MatSeqAIJSetPreallocation, MatMPIAIJSetPreallocation
    
    type(matrix_spec), intent(in) :: mat_dat
    class(ccs_matrix), allocatable, intent(out) :: M

    integer(ccs_err) :: ierr  !< Error code

    allocate(matrix_petsc :: M)

    select type (M)
      type is (matrix_petsc)

        M%modeset = .false.
        
        select type (par_env => mat_dat%par_env)
          type is(parallel_environment_mpi)

          call MatCreate(par_env%comm, M%M, ierr)

          associate(mesh => mat_dat%mesh)
            call MatSetSizes(M%M, mesh%nlocal, mesh%nlocal, PETSC_DETERMINE, PETSC_DETERMINE, ierr)
          end associate
          
          if (ierr == 0) then
            M%allocated = .true.
          end if

          call MatSetFromOptions(M%M, ierr)
          
          if (mat_dat%nnz < 1) then
            if (par_env%proc_id == par_env%root) then
              print *, "WARNING: No matrix preallocation set, potentially inefficient!"
            end if
            call MatSetUp(M%M, ierr)
          else
            call MatSeqAIJSetPreallocation(M%M, mat_dat%nnz, PETSC_NULL_INTEGER, ierr)
            call MatMPIAIJSetPreallocation(M%M, mat_dat%nnz, PETSC_NULL_INTEGER, mat_dat%nnz - 1, PETSC_NULL_INTEGER, ierr)
          end if

          class default
            print *, "Unknown parallel environment"
    
        end select

      class default
        write(*,*) "Unsupported matrix type"
        stop

    end select
    
  end subroutine

  module subroutine finalise_matrix(M)

    use petscmat, only : MAT_FINAL_ASSEMBLY
    
    class(ccs_matrix), intent(inout) :: M

    integer(ccs_err) :: ierr

    select type (M)
    type is (matrix_petsc)
       call MatAssemblyBegin(M%M, MAT_FINAL_ASSEMBLY, ierr)
       call MatAssemblyEnd(M%M, MAT_FINAL_ASSEMBLY, ierr)
    end select

  end subroutine finalise_matrix

  !> @brief Perform a parallel update of a PETSc matrix.
  !
  !> @param[in/out] M - the matrix
  module subroutine update_matrix(M)

    class(ccs_matrix), intent(inout) :: M

    select type(M)
      type is (matrix_petsc)
      
        call begin_update_matrix(M)
        call end_update_matrix(M)

      class default
        write(*,*) "Unsupported matrix type"
        stop

    end select
  
  end subroutine

  !> @brief Begin a parallel update of a PETSc matrix.
  !
  !> @details Begins the parallel update to allow overlapping comms and compute.
  !
  !> @param[in/out] M - the matrix
  module subroutine begin_update_matrix(M)

    class(ccs_matrix), intent(inout) :: M

    integer(ccs_err) :: ierr !< Error code

    select type (M)
      type is (matrix_petsc)

        call MatAssemblyBegin(M%M, MAT_FLUSH_ASSEMBLY, ierr)

      class default
        write(*,*) "Unsupported matrix type"
        stop

    end select
    
  end subroutine

  !> @brief End a parallel update of a PETSc matrix.
  !
  !> @details Ends the parallel update to allow overlapping comms and compute.
  !
  !> @param[in/out] M - the matrix
  module subroutine end_update_matrix(M)

    class(ccs_matrix), intent(inout) :: M

    integer(ccs_err) :: ierr !< Error code

    select type (M)
      type is (matrix_petsc)

        call MatAssemblyEnd(M%M, MAT_FLUSH_ASSEMBLY, ierr)

        M%modeset = .false. ! It's safe to change modes now
      class default
        write(*,*) "Unsupported matrix type"
        stop

    end select
    
  end subroutine

  module subroutine pack_one_matrix_coefficient(row_entry, col_entry, row, col, coeff, mat_coeffs)
    integer(ccs_int), intent(in) :: row_entry
    integer(ccs_int), intent(in) :: col_entry
    integer(ccs_int), intent(in) :: row
    integer(ccs_int), intent(in) :: col
    real(ccs_real), intent(in) :: coeff
    type(matrix_values), intent(inout) :: mat_coeffs

    integer(ccs_int) :: nc
    integer(ccs_int) :: validx

    mat_coeffs%row_indices(row_entry) = row - 1
    mat_coeffs%col_indices(col_entry) = col - 1

    nc = size(mat_coeffs%col_indices)

    validx = (row_entry - 1) * nc + col_entry
    mat_coeffs%values(validx) = coeff
    
  end subroutine pack_one_matrix_coefficient

  !> @brief Set values in a PETSc matrix.
  !
  !> @param[in]     mat_values - contains the values, their indices and the mode 
  !!                             to use when setting them.
  !> @param[in/out] M          - the matrix
  module subroutine set_matrix_values(mat_values, M)

    use petscmat, only : MatSetValues
    use constants, only : insert_mode, add_mode
    
    type(matrix_values), intent(in) :: mat_values
    class(ccs_matrix), intent(inout) :: M

    integer(ccs_int) :: nrows, ncols !< number of rows/columns
    integer(ccs_int) :: mode !< Add or insert values?
    
    integer(ccs_err) :: ierr !< Error code

    associate(ridx    => mat_values%row_indices, &
              cidx    => mat_values%col_indices, &
              val     => mat_values%values, &
              matmode => mat_values%setter_mode)
    
      select type (M)
        type is (matrix_petsc)

          if (M%modeset) then
            if (matmode /= M%mode) then
              print *, "ERROR: changing matrix mode without updating"
              stop 1
            end if
          else
            M%mode = matmode
            M%modeset = .true.
          end if
          
          nrows = size(ridx)
          ncols = size(cidx)
          if (nrows * ncols /= size(val)) then
            print *, "Invalid matrix values!"
            stop
          end if
          if (matmode == add_mode) then
            mode = ADD_VALUES
          else if (matmode == insert_mode) then
            mode = INSERT_VALUES
          else
            print *, "Unknown mode!"
            stop
          end if

          call MatSetValues(M%M, nrows, ridx, ncols, cidx, val, mode, ierr)

        class default
          print *, "Unknown matrix type!"
          stop

      end select

    end associate

  end subroutine

  !> @brief Set equation
  !
  !> @details Sets equations in a system of equations by zeroing out the corresponding row in the
  !!          system matrix and setting the diagonal to one such that the solution is given by
  !!          the corresponding entry in the right-hand side vector.  module subroutine set_eqn(rows, M)
  !
  !> @param[in]  rows - array of (global) row indices to set the equation on
  !> @param[in/out] M - the matrix
  module subroutine set_eqn(rows, M)

    use petsc, only : PETSC_NULL_VEC
    use petscmat, only : MatZeroRows

    integer(ccs_int), dimension(:), intent(in) :: rows
    class(ccs_matrix), intent(inout) :: M

    integer(ccs_err) :: ierr
    
    select type (M)
      type is (matrix_petsc)

        call MatZeroRows(M%M, size(rows), rows, 1.0_ccs_real, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)

      class default
        print *, "Unknown matrix type!"
        stop

    end select
    
  end subroutine

  !> @brief Perform the AXPY matrix operation using PETSc
  !
  !> @details Performs the AXPY operation
  !!          y[i] = alpha * x[i] + y[i]
  !
  !> @param[in]     alpha - a scalar value
  !> @param[in]     x     - a PETSc input matrix
  !> @param[in,out] y     - PETSc matrix serving as input, overwritten with result
  module subroutine mat_axpy(alpha, x, y)

    use petscmat, only : MatAXPY, DIFFERENT_NONZERO_PATTERN
    
    real(ccs_real), intent(in) :: alpha
    class(ccs_matrix), intent(in) :: x
    class(ccs_matrix), intent(inout) :: y

    integer(ccs_err) :: ierr !< Error code
    
    select type (x)
      type is (matrix_petsc)

        select type (y)
          type is (matrix_petsc)

            ! PETSc performs AXPY as YPAX, with result stored in Y.
            call MatAXPY(y%M, alpha, x%M, DIFFERENT_NONZERO_PATTERN, ierr)

          class default
            print *, "Unknown matrix type!"
            stop

        end select

      class default
        print *, "Unknown matrix type!"
        stop

    end select
    
  end subroutine

  !> @brief Compute the norm of a PETSc matrix
  !
  !> @param[in]  M         - the PETSc matrix
  !> @param[in]  norm_type - which norm to compute? Currently supported is the 2 norm:
  !!                         norm_type=2.
  !> @param[out] n         - the computed norm returned as the result of the function
  !!                         call.
  module function mat_norm(M, norm_type) result(n)

    use petscmat, only : NORM_1, NORM_FROBENIUS, NORM_INFINITY, MatNorm
    
    class(ccs_matrix), intent(in) :: M
    integer(ccs_int), intent(in) :: norm_type

    real(ccs_real) :: n      !< The computed norm 
    integer(ccs_err) :: ierr !< Error code
    
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
          print *, "ERROR: unknown matrix norm type ", norm_type
          stop
        end if

      class default
        print *, "Type unhandled"
        stop
    end select
    
  end function

  !> @brief Extract the diagonal elements of a matrix and store in a vector
  !
  !> @param[in]  M      - the PETSc matrix
  !> @param[out] D      - the PETSc vector containing matrix diagonal elements
  module subroutine get_matrix_diagonal(M, D)

    use petscmat, only: MatGetDiagonal

    class(ccs_matrix), intent(in)  :: M
    class(ccs_vector), intent(inout) :: D

    integer(ccs_err) :: ierr !< Error code

    select type (M)
      type is (matrix_petsc)

        select type (D)
          type is (vector_petsc)
            call MatGetDiagonal(M%M, D%v, ierr)

          class default
            print *, "Unknown vector type!"
            stop
        end select

      class default
        print *, "Unknown matrix type!"
        stop
    end select

  end subroutine

  module subroutine set_matrix_diagonal(D, M)
    use petscmat, only : MatDiagonalSet

    class(ccs_vector), intent(in) :: D
    class(ccs_matrix), intent(inout) :: M

    integer(ccs_err) :: ierr
    
    select type (M)
    type is (matrix_petsc)

      select type (D)
      type is (vector_petsc)
        call MatDiagonalSet(M%M, D%v, INSERT_VALUES, ierr)

      class default
        print *, "Unknown vector type!"
        stop
      end select

    class default
      print *, "Unknown matrix type!"
      stop
    end select

  end subroutine set_matrix_diagonal

  module subroutine zero_matrix(M)

    use petscmat, only: MatZeroEntries
    
    class(ccs_matrix), intent(inout) :: M

    integer(ccs_err) :: ierr

    select type (M)
    type is (matrix_petsc)
      call MatZeroEntries(M%M, ierr)
    class default

      print *, "Unknown matrix type!"
      stop

    end select
    
  end subroutine zero_matrix

end submodule mat_petsc
