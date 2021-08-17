!> @brief Submodule file accsvec_petsc.mod
!>
!> @details An implementation of vector objects using PETSc - the datatype and operations on it
!>          (creation, destruction, setting/getting, ...)

submodule (accsvec) accsvec_petsc

#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscvec.h>

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

    use petsc, only : PETSC_COMM_WORLD, PETSC_DECIDE
    use petscvec, only : VecCreate, VecSetSizes, VecSetFromOptions
    
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
       call VecSet(v%v, 0.0, ierr)
       v%allocated = .true.
    end select

  end subroutine

  module subroutine set_vector_values(val_dat, v)

    use petsc, only : VecSetValues
    
    use accs_types, only : vector_values
    
    class(*), intent(in) :: val_dat
    class(vector), intent(inout) :: v

    integer(accs_int) :: n
    integer(accs_err) :: ierr
    
    select type (v)
    type is (vector_petsc)
       select type (val_dat)
       type is (vector_values)
          n = size(val_dat%idx)
          call VecSetValues(v%v, n, val_dat%idx, val_dat%val, val_dat%mode, ierr)
       end select
    end select
    
  end subroutine

  module subroutine begin_update_vector(v)

    use petsc, only : VecAssemblyBegin
    
    class(vector), intent(inout) :: v

    integer(accs_err) :: ierr
    
    select type (v)
    type is (vector_petsc)
       call VecAssemblyBegin(v%v, ierr)
    end select

  end subroutine

  module subroutine end_update_vector(v)

    use petsc, only : VecAssemblyEnd
    
    class(vector), intent(inout) :: v

    integer(accs_err) :: ierr
    
    select type (v)
    type is (vector_petsc)
       call VecAssemblyEnd(v%v, ierr)
    end select

  end subroutine

  module subroutine axpy(alpha, x, y)

    use petscvec, only : VecAXPY
    
    real(accs_real), intent(in) :: alpha
    class(vector), intent(in) :: x
    class(vector), intent(inout) :: y

    integer(accs_err) :: ierr
    
    select type (x)
    type is (vector_petsc)
       select type (y)
       type is (vector_petsc)
          call VecAXPY(x%v, alpha, y%v, ierr)
       end select
    end select
    
  end subroutine

  module real(accs_real) function norm(v, norm_type)

    use petscvec, only : NORM_2, VecNorm
    
    class(vector), intent(in) :: v
    integer(accs_int), intent(in) :: norm_type

    integer(accs_err) :: ierr
    
    select type (v)
    type is (vector_petsc)
       if (norm_type == 2) then
          call VecNorm(v%v, NORM_2, norm, ierr)
       else
          print *, "ERROR: unknown vector norm type ", norm_type
          stop
       end if
    end select
    
  end function norm
  
end submodule accsvec_petsc
