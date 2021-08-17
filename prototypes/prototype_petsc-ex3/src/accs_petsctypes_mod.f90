!> @brief Module file accs_petsctypes.mod
!>
!> @details Provides petsc-extended types.

module accs_petsctypes

#include <petsc/finclude/petscvec.h>
  use petscvec, only : tVec

  use accs_types, only : vector
  
  implicit none

  private
  
  type, public, extends(vector) :: vector_petsc
     type(tVec) :: v
     logical :: allocated
   contains
     final :: free_vector_petsc
  end type vector_petsc

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
  
end module accs_petsctypes
