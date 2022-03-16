!> @brief Submodule file vec_petsc.smod
!> @build petsc
!
!> @details An implementation of vector objects using PETSc - the datatype and operations on it
!!          (creation, destruction, setting/getting, ...)
submodule (vec) vec_petsc

  use kinds, only : accs_err
  use petsctypes, only : vector_petsc
  use parallel_types_mpi, only: parallel_environment_mpi
  use constants, only: cell, face

  implicit none

contains

  !> @brief Create a PETSc-backed vector
  !
  !> @param[in]  vector_init_data vec_dat - the data describing how the vector should be created.
  !> @param[out] vector v - the vector specialised to type vector_petsc.
  module subroutine create_vector(vec_dat, v)

    use petsc, only : PETSC_DECIDE, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE, PETSC_DETERMINE
    use petscvec, only : VecCreateGhost, VecSetSizes, VecSetFromOptions, VecSet, VecSetOption, &
                         VecCreateMPI
    
    type(vector_init_data), intent(in) :: vec_dat
    class(vector), allocatable, intent(out) :: v

    integer(accs_err) :: ierr !> Error code

    allocate(vector_petsc :: v)
    
    select type (v)
      type is (vector_petsc)

      select type(par_env => vec_dat%par_env)
        type is(parallel_environment_mpi)

          associate(mesh => vec_dat%mesh)

            select case(vec_dat%loc)
              case (cell)
                call VecCreateGhost(par_env%comm, &
                   mesh%nlocal, PETSC_DECIDE, &
                   mesh%nhalo, mesh%idx_global(mesh%nlocal+1:mesh%ntotal) - 1_accs_int, &
                   v%v, ierr)
              case (face)
                call VecCreateMPI(par_env%comm, &
                   mesh%nfaces_local, PETSC_DETERMINE, v%v, ierr)
              end select
          end associate
        
          call VecSetFromOptions(v%v, ierr)
          call VecSetOption(v%v, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE, ierr)
          call VecSet(v%v, 0.0_accs_real, ierr)
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

    use petsc, only : VecSetValues, ADD_VALUES, INSERT_VALUES

    use constants, only : insert_mode, add_mode
    use types, only : vector_values
    
    class(*), intent(in) :: val_dat
    class(vector), intent(inout) :: v

    integer(accs_int) :: n    !> Number of elements to add
    integer(accs_int) :: mode !> Append or insert mode
    integer(accs_err) :: ierr !> Error code
    
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

    class(vector), intent(inout) :: v

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
    
    class(vector), intent(inout) :: v

    integer(accs_err) :: ierr !> Error code
    
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
    
    class(vector), intent(inout) :: v

    integer(accs_err) :: ierr !> Error code
    
    select type (v)
      type is (vector_petsc)

        call VecAssemblyEnd(v%v, ierr)

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

    use petsc, only : VecGhostUpdateBegin, INSERT_VALUES, SCATTER_FORWARD
    
    class(vector), intent(inout) :: v

    integer(accs_err) :: ierr !> Error code
    
    select type (v)
      type is (vector_petsc)

        call VecGhostUpdateBegin(v%v, INSERT_VALUES, SCATTER_FORWARD, ierr)

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

    use petsc, only : VecGhostUpdateEnd, INSERT_VALUES, SCATTER_FORWARD
    
    class(vector), intent(inout) :: v

    integer(accs_err) :: ierr !> Error code
    
    select type (v)
      type is (vector_petsc)

        call VecGhostUpdateEnd(v%v, INSERT_VALUES, SCATTER_FORWARD, ierr)

      class default
        print *, "Unknown vector type!"
        stop

    end select

  end subroutine

  module subroutine pack_one_vector_element(val_dat, ent, idx, val)
    type(vector_values), intent(inout) :: val_dat
    integer(accs_int), intent(in) :: ent
    integer(accs_int), intent(in) :: idx
    real(accs_real), intent(in) :: val

    val_dat%idx(ent) = idx - 1
    val_dat%val(ent) = val

  end subroutine pack_one_vector_element

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
    
    real(accs_real), intent(in) :: alpha
    class(vector), intent(in) :: x
    class(vector), intent(inout) :: y

    integer(accs_err) :: ierr !> Error code
    
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
    
    class(vector), intent(in) :: v
    integer(accs_int), intent(in) :: norm_type

    real(accs_real) :: n      !> The computed norm 
    integer(accs_err) :: ierr !> Error code
    
    n = 0.0_accs_real ! initialise norm to 0
    
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

  !> @brief Gets the data in a given vector
  !
  !> @param[in] vec   - the vector to get data from
  !> @param[in] array - an array to store the data in
  module subroutine get_vector_data(vec, array)
    use petscvec
    class(vector), intent(in) :: vec
    real(accs_real), dimension(:), pointer, intent(out) :: array
    integer :: ierr

    select type(vec)
      type is(vector_petsc)
        call VecGhostGetLocalForm(vec%v, vec%vl, ierr)
        call VecGetArrayF90(vec%vl, array, ierr)
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
    use petscvec
    class(vector), intent(in) :: vec
    real(accs_real), dimension(:), pointer, intent(in) :: array
    integer :: ierr

    select type(vec)
      type is(vector_petsc)
        call VecRestoreArrayF90(vec%vl, array, ierr)
        call VecGhostRestoreLocalForm(vec%v, vec%vl, ierr)
      class default
        print *, 'invalid vector type'
        stop
    end select
  end subroutine restore_vector_data

  module procedure zero_vector
    use petscvec
    integer(accs_err) :: ierr

    select type(vec)
    type is(vector_petsc)
      call VecZeroEntries(vec%v, ierr)
    class default
      print *, "Invalid vector type"
      stop
    end select
    
  end procedure zero_vector

  !> @brief Replaces each component of a vector by its reciprocal
  !
  !> @param[inout] vec - the vector
  module subroutine vec_reciprocal(vec)

    use petscvec, only: VecReciprocal

    class(vector), intent(inout) :: vec

    integer(accs_err) :: ierr !> Error code

    select type (vec)
      type is (vector_petsc)

       call VecReciprocal(vec%v, ierr)

      class default
        print *, "Unknown vector type!"
        stop
    end select
  end subroutine vec_reciprocal

  module procedure mult_vec_vec

    use petscvec, only : VecPointwiseMult

    integer(accs_err) :: ierr

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
      
  end procedure mult_vec_vec

  module subroutine vec_view(vec_dat, vec)

#include <petsc/finclude/petscviewer.h>

    use petscvec, only: VecView
    use petscsys
    
    type(vector_init_data), intent(in) :: vec_dat
    class(vector), intent(in) :: vec

    PetscViewer :: viewer

    integer(accs_err) :: ierr

    select type (vec)
      type is (vector_petsc)

      select type(par_env => vec_dat%par_env)
        type is(parallel_environment_mpi)

          call PetscViewerCreate(par_env%comm, viewer, ierr)
          call PetscViewerSetType(viewer, PETSCVIEWERASCII, ierr)
          call VecView(vec%v, viewer, ierr)
          call PetscViewerDestroy(viewer, ierr)

        class default
          print *, "Unknown parallel environment"
      end select

      class default
        print *, "Unknown vector type!"
        stop
    end select
  end subroutine vec_view
  
end submodule vec_petsc
