!> @brief Submodule file accsvec_petsc.mod
!>
!> @details An implementation of vector objects using PETSc - the datatype and operations on it
!>          (creation, destruction, setting/getting, ...)

submodule (accsvec) accsvec_petsc

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

    use petsc, only : PETSC_DECIDE
    use petscvec, only : VecCreate, VecSetSizes, VecSetFromOptions
    
    type(vector_init_data), intent(in) :: vec_dat
    class(vector), allocatable, intent(out) :: v
    
    integer :: ierr

    allocate(vector_petsc :: v)
    
    select type (v)
    type is (vector_petsc)
       call VecCreate(vec_dat%comm, v%v, ierr)

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

    use petsc, only : VecSetValues, ADD_VALUES, INSERT_VALUES

    use accs_constants, only : insert_mode, add_mode
    use accs_types, only : vector_values
    
    class(*), intent(in) :: val_dat
    class(vector), intent(inout) :: v

    integer(accs_int) :: n
    integer(accs_int) :: mode
    integer(accs_err) :: ierr
    
    select type (v)
    type is (vector_petsc)
       select type (val_dat)
       type is (vector_values)
          n = size(val_dat%idx)
          if (val_dat%mode == add_mode) then
             mode = ADD_VALUES
          else if (val_dat%mode == insert_mode) then
             mode = INSERT_VALUES
          else
             print *, "Unknown mode!"
             stop
          end if
          call VecSetValues(v%v, n, val_dat%idx, val_dat%val, mode, ierr)
       end select
    class default
       print *, "Unknown vector type!"
       stop
    end select
    
  end subroutine

  module subroutine update_vector(v)
    class(vector), intent(inout) :: v

    select type(v)
    type is (vector_petsc)
       call begin_update_vector(v)
       call end_update_vector(v)
    end select
  end subroutine
  
  module subroutine begin_update_vector(v)

    use petsc, only : VecAssemblyBegin
    
    class(vector), intent(inout) :: v

    integer(accs_err) :: ierr
    
    select type (v)
    type is (vector_petsc)
       call VecAssemblyBegin(v%v, ierr)
    class default
       print *, "Unknown vector type!"
       stop
    end select

  end subroutine

  module subroutine end_update_vector(v)

    use petsc, only : VecAssemblyEnd
    
    class(vector), intent(inout) :: v

    integer(accs_err) :: ierr
    
    select type (v)
    type is (vector_petsc)
       call VecAssemblyEnd(v%v, ierr)
    class default
       print *, "Unknown vector type!"
       stop
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
          ! PETSc performs AXPY as YPAX, with result stored in Y.
          call VecAXPY(y%v, alpha, x%v, ierr)
       end select
    class default
       print *, "Unknown vector type!"
       stop
    end select
    
  end subroutine

  module function norm(v, norm_type) result(n)

    use petscvec, only : NORM_2, VecNorm
    
    class(vector), intent(in) :: v
    integer(accs_int), intent(in) :: norm_type

    real(accs_real) :: n
    integer(accs_err) :: ierr
    
    n = 0.0_accs_real
    select type (v)
    type is (vector_petsc)
       if (norm_type == 2) then
          call VecNorm(v%v, NORM_2, n, ierr)
       else
          print *, "ERROR: unknown vector norm type ", norm_type
          stop
       end if
    class default
       print *, "Type unhandled"
       stop
    end select
    
  end function
  
end submodule accsvec_petsc
