submodule (mat) mat_petsc

  use kinds, only : accs_err
  use petsctypes, only : matrix_petsc
  use parallel_types_mpi, only: parallel_environment_mpi
  
  implicit none

contains

  !> @brief Create a new PETSc matrix object.
  !
  !> @param[in]  mat_dat - contains information about how the matrix should be allocated
  !> @param[out] M       - the matrix object
  module subroutine create_matrix(mat_dat, M)

    use mpi
    
    use petsc, only : PETSC_DECIDE, PETSC_NULL_INTEGER
    use petscmat, only : MatCreate, MatSetSizes, MatSetFromOptions, MatSetUp, &
                         MatSeqAIJSetPreallocation, MatMPIAIJSetPreallocation
    
    type(matrix_init_data), intent(in) :: mat_dat
    class(matrix), allocatable, intent(out) :: M

!    integer(accs_int) :: nrank !> MPI rank
    integer(accs_err) :: ierr  !> Error code

    allocate(matrix_petsc :: M)

    select type (M)
      type is (matrix_petsc)

        select type (par_env => mat_dat%par_env)
          type is(parallel_environment_mpi)

          call MatCreate(par_env%comm, M%M, ierr)
  
          if (mat_dat%rloc >= 0) then
            call MatSetSizes(M%M, mat_dat%rloc, mat_dat%cloc, PETSC_DECIDE, PETSC_DECIDE, ierr)
          else if (mat_dat%rglob >= 0) then
            call MatSetSizes(M%M, PETSC_DECIDE, PETSC_DECIDE, mat_dat%rglob, mat_dat%cglob, ierr)
          else
            print *, "ERROR: invalid matrix creation!"
            stop
          end if
          
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

    use petscmat, only : MatAssemblyBegin, MatAssemblyEnd, MAT_FINAL_ASSEMBLY
    
    class(matrix), intent(inout) :: M

    integer(accs_err) :: ierr

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

    class(matrix), intent(inout) :: M

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

    use petscmat, only : MatAssemblyBegin, MAT_FINAL_ASSEMBLY
    
    class(matrix), intent(inout) :: M

    integer(accs_err) :: ierr !> Error code

    select type (M)
      type is (matrix_petsc)

        call MatAssemblyBegin(M%M, MAT_FINAL_ASSEMBLY, ierr)

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

    use petscmat, only : MatAssemblyEnd, MAT_FINAL_ASSEMBLY
    
    class(matrix), intent(inout) :: M

    integer(accs_err) :: ierr !> Error code

    select type (M)
      type is (matrix_petsc)

        call MatAssemblyEnd(M%M, MAT_FINAL_ASSEMBLY, ierr)

      class default
        write(*,*) "Unsupported matrix type"
        stop

    end select
    
  end subroutine

  module subroutine pack_one_matrix_coefficient(mat_coeffs, row_entry, col_entry, row, col, coeff)
    type(matrix_values), intent(inout) :: mat_coeffs
    integer(accs_int), intent(in) :: row_entry
    integer(accs_int), intent(in) :: col_entry
    integer(accs_int), intent(in) :: row
    integer(accs_int), intent(in) :: col
    real(accs_real), intent(in) :: coeff

    integer(accs_int) :: nc
    integer(accs_int) :: validx

    mat_coeffs%rglob(row_entry) = row - 1
    mat_coeffs%cglob(col_entry) = col - 1

    nc = size(mat_coeffs%cglob)

    validx = (row_entry - 1) * nc + col_entry
    mat_coeffs%val(validx) = coeff
    
  end subroutine pack_one_matrix_coefficient

  !> @brief Set values in a PETSc matrix.
  !
  !> @param[in]     mat_values - contains the values, their indices and the mode 
  !!                             to use when setting them.
  !> @param[in/out] M          - the matrix
  module subroutine set_matrix_values(mat_values, M)

    use petsc, only : ADD_VALUES, INSERT_VALUES
    use petscmat, only : MatSetValues
    use constants, only : insert_mode, add_mode
    
    type(matrix_values), intent(in) :: mat_values
    class(matrix), intent(inout) :: M

    integer(accs_int) :: nrows, ncols !> number of rows/columns
    integer(accs_int) :: mode !> Add or insert values?
    
    integer(accs_err) :: ierr !> Error code

    associate(ridx    => mat_values%rglob, &
              cidx    => mat_values%cglob, &
              val     => mat_values%val, &
              matmode => mat_values%mode)
    
      select type (M)
        type is (matrix_petsc)
        
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

          call MatSetValues(M%M, nrows, ridx, ncols, cidx, real(val, kind=accs_real), mode, ierr)

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

    integer(accs_int), dimension(:), intent(in) :: rows
    class(matrix), intent(inout) :: M

    integer(accs_err) :: ierr
    
    select type (M)
      type is (matrix_petsc)

        call MatZeroRows(M%M, size(rows), rows, 1.0_accs_real, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)

      class default
        print *, "Unknown matrix type!"
        stop

    end select
    
  end subroutine

end submodule mat_petsc
