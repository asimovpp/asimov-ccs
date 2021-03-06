!>  Submodule file vec_petsc.smod
!
!>  @build petsc
!
!v  An implementation of vector objects using PETSc - the datatype and operations on it
!   (creation, destruction, setting/getting, ...)
submodule (vec) vec_petsc
#include "ccs_macros.inc"

  use kinds, only : ccs_err
  use petsctypes, only : vector_petsc
  use parallel_types_mpi, only: parallel_environment_mpi
  use constants, only: cell, face
  use petsc, only: ADD_VALUES, INSERT_VALUES, SCATTER_FORWARD
  use utils, only: debug_print

  implicit none

contains

  !>  Create a PETSc-backed vector
  module subroutine create_vector(vec_properties, v)

    use petsc, only : PETSC_DECIDE, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE
    use petscvec, only : VecCreateGhost, VecSetSizes, VecSetFromOptions, VecSet, VecSetOption, &
                         VecCreate
    
    type(vector_spec), intent(in) :: vec_properties     !< the data describing how the vector should be created.
    class(ccs_vector), allocatable, intent(out) :: v    !< the vector specialised to type vector_petsc.

    integer(ccs_err) :: ierr !< Error code

    allocate(vector_petsc :: v)
    
    select type (v)
      type is (vector_petsc)

      v%modeset = .false.
        
      select type(par_env => vec_properties%par_env)
        type is(parallel_environment_mpi)

          associate(mesh => vec_properties%mesh)

            select case(vec_properties%storage_location)
            case (cell)
              call VecCreateGhost(par_env%comm, &
                   mesh%nlocal, PETSC_DECIDE, &
                   mesh%nhalo, &
                   mesh%global_indices(min(mesh%nlocal+1, mesh%ntotal):mesh%ntotal) - 1_ccs_int, &
                   v%v, ierr)
              ! Vector has ghost points, store this information
              v%ghosted = .true.
            case (face)
              call VecCreate(par_env%comm, v%v, ierr)
              call VecSetSizes(v%v, mesh%nfaces_local, PETSC_DECIDE, ierr)

              ! Vector doesn't have ghost points, store this information
              v%ghosted = .false.
            end select
          end associate
        
          call VecSetFromOptions(v%v, ierr)
          call VecSetOption(v%v, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE, ierr)
          call VecSet(v%v, 0.0_ccs_real, ierr)
          v%allocated = .true.

        class default
          print *, "Unknown parallel environment"
    
        end select

      class default
        print *, "Unknown vector type!"
        stop

    end select

  end subroutine

  !>  Sets values in a PETSc vector
  module subroutine set_vector_values(val_dat, v)

    use petsc, only : VecSetValues

    use constants, only : insert_mode, add_mode
    
    class(*), intent(in) :: val_dat         !< contains the values, their indices and the mode to use for setting them.
    class(ccs_vector), intent(inout) :: v   !< the PETSc vector.  
                                            
    integer(ccs_int) :: n    !< Number of elements to add
    integer(ccs_int) :: mode !< Append or insert mode
    integer(ccs_err) :: ierr !< Error code
    
    select type (v)
      type is (vector_petsc)
      
        select type (val_dat)
          type is (vector_values)

            ! First check if safe to set
            if (v%modeset) then
              if (val_dat%setter_mode /= v%mode) then
                print *, "ERROR: trying to set vector using different mode without updating!"
                stop 1
              end if
            else
              v%mode = val_dat%setter_mode
              v%modeset = .true.
            end if
              
            n = size(val_dat%global_indices)
            if (val_dat%setter_mode == add_mode) then
              mode = ADD_VALUES
            else if (val_dat%setter_mode == insert_mode) then
              mode = INSERT_VALUES
            else
              print *, "Unknown mode!"
              stop
            end if

            call VecSetValues(v%v, n, val_dat%global_indices, val_dat%values, mode, ierr)

          class default
            print *, "Unknown vector value type!"
            stop
   
        end select

      class default
        print *, "Unknown vector type!"
        stop

    end select
    
  end subroutine

  !>  Perform a parallel update of a PETSc vector
  module subroutine update_vector(v)

    class(ccs_vector), intent(inout) :: v   !< the PETSc vector

    select type(v)
      type is (vector_petsc)

        call begin_update_vector(v)
        call end_update_vector(v)

        call begin_ghost_update_vector(v)
        call end_ghost_update_vector(v)
        
      class default
        print *, "Unknown vector type!"
        stop

    end select

  end subroutine

  !>  Begin a parallel update of a PETSc vector
  !
  !>  Begins the parallel update to allow overlapping comms and compute
  module subroutine begin_update_vector(v)

    use petsc, only : VecAssemblyBegin
    
    class(ccs_vector), intent(inout) :: v   !< the PETSc vector

    integer(ccs_err) :: ierr !< Error code
    
    select type (v)
      type is (vector_petsc)

        call VecAssemblyBegin(v%v, ierr)

      class default
        print *, "Unknown vector type!"
        stop

    end select

  end subroutine

  !>  End a parallel update of a PETSc vector.
  !
  !>  Ends the parallel update to allow overlapping comms and compute.
  module subroutine end_update_vector(v)

    use petsc, only : VecAssemblyEnd
    
    class(ccs_vector), intent(inout) :: v   !< the PETSc vector

    integer(ccs_err) :: ierr !< Error code
    
    select type (v)
      type is (vector_petsc)

        call VecAssemblyEnd(v%v, ierr)

        v%modeset = .false. ! It is now safe to change value setting mode
        
      class default
        print *, "Unknown vector type!"
        stop

    end select

  end subroutine

  !>  Begin a ghost update of a PETSc vector
  !
  !>  Begins the ghost update to allow overlapping comms and compute
  subroutine begin_ghost_update_vector(v)

    use petsc, only : VecGhostUpdateBegin
    
    class(ccs_vector), intent(inout) :: v   !< the PETSc vector

    integer(ccs_err) :: ierr !< Error code
    
    select type (v)
      type is (vector_petsc)

        if (v%ghosted) then
          ! Cant update ghosts if not ghost points!
          call VecGhostUpdateBegin(v%v, INSERT_VALUES, SCATTER_FORWARD, ierr)
        end if
        
      class default
        print *, "Unknown vector type!"
        stop

    end select

  end subroutine

  !>  End a ghost update of a PETSc vector.
  !
  !>  Ends the ghost update to allow overlapping comms and compute.
  subroutine end_ghost_update_vector(v)

    use petsc, only : VecGhostUpdateEnd
    
    class(ccs_vector), intent(inout) :: v   !< the PETSc vector

    integer(ccs_err) :: ierr !< Error code
    
    select type (v)
      type is (vector_petsc)

        if (v%ghosted) then
          ! Cant update ghosts if not ghost points!
          call VecGhostUpdateEnd(v%v, INSERT_VALUES, SCATTER_FORWARD, ierr)
        end if
        
      class default
        print *, "Unknown vector type!"
        stop

    end select

  end subroutine

  module subroutine pack_one_vector_element(ent, idx, val, val_dat)
    integer(ccs_int), intent(in) :: ent
    integer(ccs_int), intent(in) :: idx
    real(ccs_real), intent(in) :: val
    type(vector_values), intent(inout) :: val_dat

    val_dat%global_indices(ent) = idx - 1
    val_dat%values(ent) = val

  end subroutine pack_one_vector_element

  !>  Perform the AXPY vector operation using PETSc
  !
  !v  Performs the AXPY operation
  !          y[i] = alpha * x[i] + y[i]
  module subroutine vec_axpy(alpha, x, y)

    use petscvec, only : VecAXPY
    
    real(ccs_real), intent(in) :: alpha     !< a scalar value
    class(ccs_vector), intent(in) :: x      !< a PETSc input vector
    class(ccs_vector), intent(inout) :: y   !< PETSc vector serving as input, overwritten with result

    integer(ccs_err) :: ierr !< Error code
    
    select type (x)
      type is (vector_petsc)

        select type (y)
          type is (vector_petsc)

            ! PETSc performs AXPY as YPAX, with result stored in Y.
            call VecAXPY(y%v, alpha, x%v, ierr)

          class default
            print *, "Unknown vector type!"
            stop

        end select

      class default
        print *, "Unknown vector type!"
        stop

    end select
    
  end subroutine

  !>  Compute the norm of a PETSc vector
  module function vec_norm(v, norm_type) result(n)

    use petscvec, only : NORM_2, VecNorm
    
    class(ccs_vector), intent(in) :: v          !< the PETSc vector
    integer(ccs_int), intent(in) :: norm_type   !< which norm to compute? Currently supported is the 2 norm: norm_type=2.

    real(ccs_real) :: n      !< The computed norm 
    integer(ccs_err) :: ierr !< Error code
    
    n = 0.0_ccs_real ! initialise norm to 0
    
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

  !>  Gets the data in a given vector
  !
  !> @param[in] vec   - the vector to get data from
  !> @param[in] array - an array to store the data in
  module subroutine get_vector_data(vec, array)
    use petscvec, only: VecGhostGetLocalForm, VecGetArrayF90
    class(ccs_vector), intent(in) :: vec
    real(ccs_real), dimension(:), pointer, intent(out) :: array
    integer :: ierr

    select type(vec)
    type is(vector_petsc)
      if (vec%modeset) then
        print *, "WARNING: trying to access vector without updating"
      end if
      
      if (vec%ghosted) then
        call VecGhostGetLocalForm(vec%v, vec%v_local, ierr)
        call VecGetArrayF90(vec%v_local, array, ierr)
      else
        call VecGetArrayF90(vec%v, array, ierr)
      end if
    class default
      print *, 'invalid vector type'
      stop
    end select
  end subroutine get_vector_data

  !>  Resets the vector data if required for further processing
  !
  !> @param[in] vec   - the vector to reset
  !> @param[in] array - the array containing the data to restore
  module subroutine restore_vector_data(vec, array)
    use petscvec, only: VecRestoreArrayF90, VecGhostRestoreLocalForm

    class(ccs_vector), intent(in) :: vec
    real(ccs_real), dimension(:), pointer, intent(in) :: array

    integer :: ierr

    select type(vec)
    type is(vector_petsc)
      if (vec%ghosted) then
        call VecRestoreArrayF90(vec%v_local, array, ierr)
        call VecGhostRestoreLocalForm(vec%v, vec%v_local, ierr)
      else
        call VecRestoreArrayF90(vec%v, array, ierr)
      end if
    class default
      print *, 'invalid vector type'
      stop
    end select
  end subroutine restore_vector_data

  module subroutine zero_vector(vec)
    use petscvec, only: VecZeroEntries

    class(ccs_vector), intent(inout) :: vec

    integer(ccs_err) :: ierr

    select type(vec)
    type is(vector_petsc)
      call VecZeroEntries(vec%v, ierr)
    class default
      print *, "Invalid vector type"
      stop
    end select
    
  end subroutine zero_vector

  !>  Replaces each component of a vector by its reciprocal
  !
  !> @param[inout] vec - the vector
  module subroutine vec_reciprocal(vec)

    use petscvec, only: VecReciprocal

    class(ccs_vector), intent(inout) :: vec

    integer(ccs_err) :: ierr !< Error code

    select type (vec)
      type is (vector_petsc)

       call VecReciprocal(vec%v, ierr)

      class default
        print *, "Unknown vector type!"
        stop
    end select
  end subroutine vec_reciprocal

  module subroutine mult_vec_vec(a, b)
    use petscvec, only : VecPointwiseMult

    class(ccs_vector), intent(in) :: a
    class(ccs_vector), intent(inout) :: b

    integer(ccs_err) :: ierr

    select type(a)
    type is (vector_petsc)
      select type(b)
      type is (vector_petsc)
        call VecPointwiseMult(b%v, a%v, b%v, ierr)
      class default
        print *, "Unknown vector type!"
        stop
      end select
    class default
      print *, "Unknown vector type!"
      stop
    end select
      
  end subroutine mult_vec_vec

  module subroutine scale_vec(alpha, v)
    use petscvec, only : VecScale

    real(ccs_real), intent(in) :: alpha
    class(ccs_vector), intent(inout) :: v
  
    integer(ccs_err) :: ierr

    select type(v)
    type is(vector_petsc)
      call VecScale(v%v, alpha, ierr)
    class default
      print *, "Unknown vector type!"
      stop
    end select
  
  end subroutine

!   module subroutine vec_view(vec_properties, vec)

! #include <petsc/finclude/petscviewer.h>

!     use petscvec, only: VecView
!     use petscsys
    
!     type(vector_spec), intent(in) :: vec_properties
!     class(ccs_vector), intent(in) :: vec

!     PetscViewer :: viewer

!     integer(ccs_err) :: ierr

!     select type (vec)
!       type is (vector_petsc)

!       select type(par_env => vec_properties%par_env)
!         type is(parallel_environment_mpi)

!           call PetscViewerCreate(par_env%comm, viewer, ierr)
!           call PetscViewerSetType(viewer, PETSCVIEWERASCII, ierr)
!           call VecView(vec%v, viewer, ierr)
!           call PetscViewerDestroy(viewer, ierr)

!         class default
!           print *, "Unknown parallel environment"
!       end select

!       class default
!         print *, "Unknown vector type!"
!         stop
!     end select
!   end subroutine vec_view

end submodule vec_petsc
