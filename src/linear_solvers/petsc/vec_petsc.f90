!> @brief Submodule file vec_petsc.smod
!> @build petsc
!
!> @details An implementation of vector objects using PETSc - the datatype and operations on it
!!          (creation, destruction, setting/getting, ...)
submodule (vec) vec_petsc

  use kinds, only : ccs_err
  use petsctypes, only : vector_petsc
  use parallel_types_mpi, only: parallel_environment_mpi
  use constants, only: cell, face
  use petsc, only: ADD_VALUES, INSERT_VALUES, SCATTER_FORWARD

  implicit none

contains

  !> @brief Create a PETSc-backed vector
  !
  !> @param[in]  vector_spec vec_dat - the data describing how the vector should be created.
  !> @param[out] vector v - the vector specialised to type vector_petsc.
  module subroutine create_vector(vec_dat, v)

    use petsc, only : PETSC_DECIDE, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE
    use petscvec, only : VecCreateGhost, VecSetSizes, VecSetFromOptions, VecSet, VecSetOption, &
                         VecCreate
    
    type(vector_spec), intent(in) :: vec_dat
    class(ccs_vector), allocatable, intent(out) :: v

    integer(ccs_int), dimension(:), allocatable :: global_halo_indices
    integer(ccs_int) :: i
    integer(ccs_err) :: ierr !> Error code

    allocate(vector_petsc :: v)
    
    select type (v)
      type is (vector_petsc)

      v%modeset = .false.
        
      select type(par_env => vec_dat%par_env)
        type is(parallel_environment_mpi)

          associate(mesh => vec_dat%mesh)

            select case(vec_dat%storage_location)
            case (cell)
              associate(nhalo => mesh%nhalo, &
                   nlocal => mesh%nlocal, &
                   idx_global => mesh%idx_global)
                allocate(global_halo_indices(nhalo))
                do i = 1, nhalo
                  global_halo_indices(i) = idx_global(i + nlocal) - 1_ccs_int
                end do
                call VecCreateGhost(par_env%comm, &
                     nlocal, PETSC_DECIDE, &
                     nhalo, global_halo_indices, &
                     v%v, ierr)
                deallocate(global_halo_indices)
              end associate
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

  !> @brief Sets values in a PETSc vector
  !
  !> @param[in] val_dat - contains the values, their indices and the mode to use for setting
  !!                      them.
  !> @param[in,out] v   - the PETSc vector.
  module subroutine set_vector_values(val_dat, v)

    use petsc, only : VecSetValues

    use constants, only : insert_mode, add_mode
    
    class(*), intent(in) :: val_dat
    class(ccs_vector), intent(inout) :: v

    integer(ccs_int) :: n    !> Number of elements to add
    integer(ccs_int) :: mode !> Append or insert mode
    integer(ccs_err) :: ierr !> Error code
    
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
              
            n = size(val_dat%indices)
            if (val_dat%setter_mode == add_mode) then
              mode = ADD_VALUES
            else if (val_dat%setter_mode == insert_mode) then
              mode = INSERT_VALUES
            else
              print *, "Unknown mode!"
              stop
            end if

            call VecSetValues(v%v, n, val_dat%indices, val_dat%values, mode, ierr)

          class default
            print *, "Unknown vector value type!"
            stop
   
        end select

      class default
        print *, "Unknown vector type!"
        stop

    end select
    
  end subroutine

  !> @brief Perform a parallel update of a PETSc vector
  !
  !> @param[in,out] v - the PETSc vector
  module subroutine update_vector(v)

    class(ccs_vector), intent(inout) :: v

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

  !> @brief Begin a parallel update of a PETSc vector
  !
  !> @details Begins the parallel update to allow overlapping comms and compute
  !
  !> @param[in,out] v - the PETSc vector
  module subroutine begin_update_vector(v)

    use petsc, only : VecAssemblyBegin
    
    class(ccs_vector), intent(inout) :: v

    integer(ccs_err) :: ierr !> Error code
    
    select type (v)
      type is (vector_petsc)

        call VecAssemblyBegin(v%v, ierr)

      class default
        print *, "Unknown vector type!"
        stop

    end select

  end subroutine

  !> @brief End a parallel update of a PETSc vector.
  !
  !> @details Ends the parallel update to allow overlapping comms and compute.
  !
  !> @param[in,out] v - the PETSc vector
  module subroutine end_update_vector(v)

    use petsc, only : VecAssemblyEnd
    
    class(ccs_vector), intent(inout) :: v

    integer(ccs_err) :: ierr !> Error code
    
    select type (v)
      type is (vector_petsc)

        call VecAssemblyEnd(v%v, ierr)

        v%modeset = .false. ! It is now safe to change value setting mode
        
      class default
        print *, "Unknown vector type!"
        stop

    end select

  end subroutine

  !> @brief Begin a ghost update of a PETSc vector
  !
  !> @details Begins the ghost update to allow overlapping comms and compute
  !
  !> @param[in,out] v - the PETSc vector
  subroutine begin_ghost_update_vector(v)

    use petsc, only : VecGhostUpdateBegin
    
    class(ccs_vector), intent(inout) :: v

    integer(ccs_err) :: ierr !> Error code
    
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

  !> @brief End a ghost update of a PETSc vector.
  !
  !> @details Ends the ghost update to allow overlapping comms and compute.
  !
  !> @param[in,out] v - the PETSc vector
  subroutine end_ghost_update_vector(v)

    use petsc, only : VecGhostUpdateEnd
    
    class(ccs_vector), intent(inout) :: v

    integer(ccs_err) :: ierr !> Error code
    
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

  !> @brief Perform the AXPY vector operation using PETSc
  !
  !> @details Performs the AXPY operation
  !!          y[i] = alpha * x[i] + y[i]
  !
  !> @param[in]     alpha - a scalar value
  !> @param[in]     x     - a PETSc input vector
  !> @param[in,out] y     - PETSc vector serving as input, overwritten with result
  module subroutine vec_axpy(alpha, x, y)

    use petscvec, only : VecAXPY
    
    real(ccs_real), intent(in) :: alpha
    class(ccs_vector), intent(in) :: x
    class(ccs_vector), intent(inout) :: y

    integer(ccs_err) :: ierr !> Error code
    
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

  !> @brief Compute the norm of a PETSc vector
  !
  !> @param[in]  v         - the PETSc vector
  !> @param[in]  norm_type - which norm to compute? Currently supported is the 2 norm:
  !!                         norm_type=2.
  !> @param[out] n         - the computed norm returned as the result of the function
  !!                         call.
  module function vec_norm(v, norm_type) result(n)

    use petscvec, only : NORM_2, VecNorm
    
    class(ccs_vector), intent(in) :: v
    integer(ccs_int), intent(in) :: norm_type

    real(ccs_real) :: n      !> The computed norm 
    integer(ccs_err) :: ierr !> Error code
    
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

  module procedure clear_vector_values_entries

    val_dat%indices(:) = -1 ! PETSc ignores -ve indices, used as "empty" indicator
    val_dat%values(:) = 0.0_ccs_real
    
  end procedure clear_vector_values_entries
  
  module procedure set_vector_values_row

    integer(ccs_int), dimension(rank(val_dat%indices)) :: idxs
    integer(ccs_int) :: i
    integer(ccs_int) :: petsc_row

    petsc_row = row - 1 ! PETSc is zero-indexed
    
    idxs = findloc(val_dat%indices, petsc_row, kind=ccs_int)
    i = idxs(1) ! We want the first entry
    if (i == 0) then
      ! New entry
      idxs = findloc(val_dat%indices, -1_ccs_int, kind=ccs_int)
      i = idxs(1) ! We want the first entry
      if (i == 0) then
        print *, "ERROR: Couldn't find a free entry in vector values!"
        stop
      end if
    end if
    
    val_dat%current_entry = i
    val_dat%indices(i) = petsc_row
    
  end procedure set_vector_values_row

  !> @brief Gets the data in a given vector
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
        call VecGhostGetLocalForm(vec%v, vec%vl, ierr)
        call VecGetArrayF90(vec%vl, array, ierr)
      else
        call VecGetArrayF90(vec%v, array, ierr)
      end if
    class default
      print *, 'invalid vector type'
      stop
    end select
  end subroutine get_vector_data

  !> @brief Resets the vector data if required for further processing
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
        call VecRestoreArrayF90(vec%vl, array, ierr)
        call VecGhostRestoreLocalForm(vec%v, vec%vl, ierr)
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

  !> @brief Replaces each component of a vector by its reciprocal
  !
  !> @param[inout] vec - the vector
  module subroutine vec_reciprocal(vec)

    use petscvec, only: VecReciprocal

    class(ccs_vector), intent(inout) :: vec

    integer(ccs_err) :: ierr !> Error code

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

!   module subroutine vec_view(vec_dat, vec)

! #include <petsc/finclude/petscviewer.h>

!     use petscvec, only: VecView
!     use petscsys
    
!     type(vector_spec), intent(in) :: vec_dat
!     class(ccs_vector), intent(in) :: vec

!     PetscViewer :: viewer

!     integer(ccs_err) :: ierr

!     select type (vec)
!       type is (vector_petsc)

!       select type(par_env => vec_dat%par_env)
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
