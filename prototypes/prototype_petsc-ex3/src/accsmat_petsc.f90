submodule (accsmat) accsmat_petsc

  use accs_types, only : matrix_init_data
  use accs_petsctypes, only : matrix_petsc
  
  implicit none

contains

  module subroutine create_matrix(mat_dat, M)

    use petsc, only : PETSC_DECIDE
    use petscmat, only : MatCreate, MatSetSizes, MatSetFromOptions, MatSetUp
    
    type(matrix_init_data), intent(in) :: mat_dat
    class(matrix), allocatable, intent(out) :: M

    integer :: ierr

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

       call MatSetFromOptions(M%M, ierr)
       call MatSetUp(M%M, ierr)
    end select
    
  end subroutine

  module subroutine begin_update_matrix(M)

    use petscmat, only : MatAssemblyBegin, MAT_FINAL_ASSEMBLY
    
    class(matrix), intent(inout) :: M

    integer :: ierr

    select type (M)
    type is (matrix_petsc)
       call MatAssemblyBegin(M%M, MAT_FINAL_ASSEMBLY, ierr)
    end select
    
  end subroutine

  module subroutine end_update_matrix(M)

    use petscmat, only : MatAssemblyEnd, MAT_FINAL_ASSEMBLY
    
    class(matrix), intent(inout) :: M

    integer :: ierr

    select type (M)
    type is (matrix_petsc)
       call MatAssemblyEnd(M%M, MAT_FINAL_ASSEMBLY, ierr)
    end select
    
  end subroutine

  module subroutine set_matrix_values(mat_values, M)
    type(matrix_values), intent(in) :: mat_values
    class(matrix), intent(inout) :: M

  end subroutine

end submodule accsmat_petsc
