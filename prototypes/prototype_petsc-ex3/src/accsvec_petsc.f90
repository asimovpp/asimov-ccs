!> @brief Submodule file accsvec_petsc.mod
!>
!> @details An implementation of vector objects using PETSc - the datatype and operations on it
!>          (creation, destruction, setting/getting, ...)

submodule (accsvec) accsvec_petsc

#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscvec.h>
  use petsc, only : PETSC_COMM_WORLD, PETSC_DECIDE
  use petscvec, only : VecCreate, VecSetSizes, VecSetFromOptions

  use accs_kinds, only : accs_int, accs_err
  use accs_types, only : vector, vector_init_data
  use accs_petsctypes, only : vector_petsc

  implicit none

contains

  module subroutine create_vector(vec_dat, v)
    !> @brief Creates a PETSc-backed vector.
    !>
    !> @param[in] vector_innit_data vec_dat - the data describing how the vector should be created.
    !> @param[out] vector v - the vector specialised to type vector_petsc.
    
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
       v%allocated = .true.
    end select

  end subroutine

  module subroutine free_vector(v)
    !> @brief Destroys a PETSc-backed vector.
    !>
    !> @param[in] vector v - the vector to be destroyed.
    
    class(vector), allocatable, intent(inout) :: v

    integer :: ierr

    if (allocated(v)) then
       select type (v)
       type is (vector_petsc)
          
          if (v%allocated) then
             call VecDestroy(v%v, ierr)
             v%allocated = .false.
          else
             print *, "WARNING: attempted double free of vector"
          end if

          ! deallocate(v) ! XXX: I feel like we should deallocate(v) here, but it won't compile...
       end select
    else
       print *, "WARNING: attempted double free of vector"
    end if
    
  end subroutine

  module subroutine set_vector_values(val_dat, v)

    use accs_types, only : vector_values
    
    class(*), intent(in) :: val_dat
    class(vector), intent(inout) :: v

    integer(accs_int) :: n
    
    select type (v)
    type is (vector_petsc)
       select type (val_dat)
       type is (vector_values)
          n = size(val_dat%idx)
          call VecSetValues(v%v, n, val_dat%idx, val_dat%val, val_dat%val)
       end select
    end select
    
  end subroutine

  module subroutine begin_update_vector(v)

    class(vector), intent(inout) :: v

    integer(accs_err) :: ierr
    
    select type (v)
    type is (vector_petsc)
       call VecAssemblyBegin(v, ierr)
    end select

  end subroutine

  module subroutine end_update_vector(v)

    class(vector), intent(inout) :: v

    integer(accs_err) :: ierr
    
    select type (v)
    type is (vector_petsc)
       call VecAssemblyEnd(v, ierr)
    end select

  end subroutine
  
end submodule accsvec_petsc
