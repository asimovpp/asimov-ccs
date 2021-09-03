submodule (accsmat) accsmat_petsc

  use accs_kinds, only : accs_int, accs_real, accs_err
  use accs_types, only : matrix, matrix_init_data
  use accs_petsctypes, only : matrix_petsc
  
  implicit none

contains

  module subroutine create_matrix(mat_dat, M)

    use petsc, only : PETSC_DECIDE
    use petscmat, only : MatCreate, MatSetSizes, MatSetFromOptions, MatSetUp
    
    type(matrix_init_data), intent(in) :: mat_dat
    class(matrix), allocatable, intent(out) :: M

    integer(accs_err) :: ierr

    allocate(matrix_petsc :: M)
    select type (M)
    type is (matrix_petsc)
       call MatCreate(mat_dat%comm, M%M, ierr)
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
       call MatSetUp(M%M, ierr)
    end select
    
  end subroutine

  module subroutine update_matrix(M)
    class(matrix), intent(inout) :: M

    select type(M)
    type is (matrix_petsc)
       call begin_update_matrix(M)
       call end_update_matrix(M)
    end select
  end subroutine
  
  module subroutine begin_update_matrix(M)

    use petscmat, only : MatAssemblyBegin, MAT_FINAL_ASSEMBLY
    
    class(matrix), intent(inout) :: M

    integer(accs_err) :: ierr

    select type (M)
    type is (matrix_petsc)
       call MatAssemblyBegin(M%M, MAT_FINAL_ASSEMBLY, ierr)
    end select
    
  end subroutine

  module subroutine end_update_matrix(M)

    use petscmat, only : MatAssemblyEnd, MAT_FINAL_ASSEMBLY
    
    class(matrix), intent(inout) :: M

    integer(accs_err) :: ierr

    select type (M)
    type is (matrix_petsc)
       call MatAssemblyEnd(M%M, MAT_FINAL_ASSEMBLY, ierr)
    end select
    
  end subroutine

  module subroutine set_matrix_values(mat_values, M)

    use petsc, only : ADD_VALUES, INSERT_VALUES
    use petscmat, only : MatSetValues
    use accs_constants, only : insert_mode, add_mode
    
    type(matrix_values), intent(in) :: mat_values
    class(matrix), intent(inout) :: M

    integer(accs_int) :: nrows, ncols
    integer(accs_int) :: mode
    
    integer(accs_err) :: ierr

    associate(ridx=>mat_values%rglob, &
         cidx=>mat_values%cglob, &
         val=>mat_values%val, &
         matmode=>mat_values%mode)
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
         call MatSetValues(M%M, nrows, ridx, ncols, cidx, val, mode, ierr)
      class default
         print *, "Unknown matrix type!"
         stop
      end select
    end associate

  end subroutine

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

end submodule accsmat_petsc