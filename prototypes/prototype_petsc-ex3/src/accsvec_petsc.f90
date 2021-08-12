!> @brief Submodule file accsvec_petsc.mod
!>
!> @details An implementation of vector objects using PETSc - the datatype and operations on it
!>          (creation, destruction, setting/getting, ...)

submodule (accsvec) accsvec_petsc

#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscvec.h>
  use petsc, only : PETSC_COMM_WORLD, PETSC_DECIDE
  use petscvec, only : tVec, VecCreate, VecSetSizes, VecSetFromOptions

  implicit none
  
  type, extends(vector) :: vector_petsc
     type(tVec) :: v
  end type vector_petsc

contains

  module subroutine create_vector(vec_dat, v)

    type(vector_init_data), intent(in) :: vec_dat
    class(vector), allocatable, intent(out) :: v
    
    integer :: ierr

    allocate(vector_petsc :: v)
    
    select type (v)
    type is (vector_petsc)
       call VecCreate(PETSC_COMM_WORLD, v%v, ierr)

       if (vec_dat%nloc >= 0) then
          call VecSetSizes(v%v, vec_dat%nloc, PETSC_DECIDE, ierr)
       else if (vec_dat%nglob > 0) then
          call VecSetSizes(v%v, PETSC_DECIDE, vec_dat%nglob, ierr)
       else
          print *, "ERROR: invalid vector creation!"
          stop
       end if
       
       call VecSetFromOptions(v%v, ierr)
    end select

  end subroutine
  
end submodule accsvec_petsc
