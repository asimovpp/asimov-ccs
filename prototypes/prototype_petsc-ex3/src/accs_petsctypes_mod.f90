!> @brief Module file accs_petsctypes.mod
!>
!> @details Provides petsc-extended types.

module accs_petsctypes

  use petscvec, only : tVec
  use petscmat, only : tMat

  use accs_types, only : vector, matrix
  
  implicit none

  private
  
  type, public, extends(vector) :: vector_petsc
     type(tVec) :: v
     logical :: allocated
   contains
     final :: free_vector_petsc
  end type vector_petsc

  type, public, extends(matrix) :: matrix_petsc
     type(tMat) :: M
     logical :: allocated
   contains
     final :: free_matrix_petsc
  end type matrix_petsc

contains
  
  module subroutine free_vector_petsc(v)
    !> @brief Destroys a PETSc-backed vector.
    !>
    !> @param[in] vector v - the vector to be destroyed.
    
    type(vector_petsc), intent(inout) :: v

    integer :: ierr

    if (v%allocated) then
       call VecDestroy(v%v, ierr)
       v%allocated = .false.
    else
       print *, "WARNING: attempted double free of vector"
    end if
    
  end subroutine

  module subroutine free_matrix_petsc(M)
    !> @brief Destroys a PETSc-backed matrix.
    !>
    !> @param[in] matrix M - the matrix to be destroyed.

    type(matrix_petsc), intent(inout) :: M

    integer :: ierr

    if (M%allocated) then
       call MatDestroy(M%M, ierr)
       M%allocated = .false.
    else
       print *, "WARNING: attempted double free of matrix"
    end if

  end subroutine
  
end module accs_petsctypes
