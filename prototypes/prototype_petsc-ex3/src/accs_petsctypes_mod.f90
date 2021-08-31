!> @brief Module file accs_petsctypes.mod
!>
!> @details Provides petsc-extended types.

module accs_petsctypes

  use petscvec, only : tVec
  use petscmat, only : tMat

  use accs_types, only : vector, matrix
  
  implicit none

  private

  !> @brief Implements the vector class backed by a PETSc vector
  type, public, extends(vector) :: vector_petsc
     type(tVec) :: v      !> The PETSc vector
     logical :: allocated !> Indicates whether the PETSc vector has been allocated
   contains
     final :: free_vector_petsc
  end type vector_petsc

  !> @brief Implements the matrix class backed by a PETSc matrix
  type, public, extends(matrix) :: matrix_petsc
     type(tMat) :: M      !> The PETSc matrix
     logical :: allocated !> Indicates whether the PETSc matrix has been allocated
   contains
     final :: free_matrix_petsc
  end type matrix_petsc

contains
  
  !> @brief Destroys a PETSc-backed vector.
  !>
  !> @param[in] vector v - the vector to be destroyed.
  !>
  !> @details Destructor called by deallocating a vector_petsc - confirms the PETSc vector object is
  !>          allocated and calls the necessary destructor on the wrapped PETSc vector object, sets
  !>          the allocated flag to .false. to prevent double free's.
  module subroutine free_vector_petsc(v)
    
    type(vector_petsc), intent(inout) :: v

    integer :: ierr

    if (v%allocated) then
       call VecDestroy(v%v, ierr)
       v%allocated = .false.
    else
       print *, "WARNING: attempted double free of vector"
    end if
    
  end subroutine
  
  !> @brief Destroys a PETSc-backed matrix.
  !>
  !> @param[in] matrix M - the matrix to be destroyed.
  !>
  !> @details Destructor called by deallocating a vector_petsc - confirms the PETSc vector object is
  !>          allocated and calls the necessary destructor on the wrapped PETSc vector object, sets
  !>          the allocated flag to .false. to prevent double free's.
  module subroutine free_matrix_petsc(M)

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
