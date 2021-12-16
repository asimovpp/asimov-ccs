!> @brief Submodule file vec_petsc.smod
!> @build petsc
!
!> @details An implementation of vector objects using PETSc - the datatype and operations on it
!!          (creation, destruction, setting/getting, ...)
submodule (vec) vec_petsc

  use kinds, only : accs_err
  use petsctypes, only : vector_petsc
  use parallel_types_mpi, only: parallel_environment_mpi

  implicit none

contains

  !> @brief Create a PETSc-backed vector
  !
  !> @param[in]  vector_innit_data vec_dat - the data describing how the vector should be created.
  !> @param[out] vector v - the vector specialised to type vector_petsc.
  module subroutine create_vector(vec_dat, v)

    use petsc, only : PETSC_DECIDE, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE
    use petscvec, only : VecCreate, VecSetSizes, VecSetFromOptions, VecSet, VecSetOption
    
    type(vector_init_data), intent(in) :: vec_dat
    class(vector), allocatable, intent(out) :: v

    integer(accs_err) :: ierr !> Error code

    allocate(vector_petsc :: v)
    
    select type (v)
      type is (vector_petsc)

        select type(par_env => vec_dat%par_env)
        type is(parallel_environment_mpi)
          call VecCreate(par_env%comm, v%v, ierr)

          if (vec_dat%nloc >= 0) then
            call VecSetSizes(v%v, vec_dat%nloc, PETSC_DECIDE, ierr)
          else if (vec_dat%nglob > 0) then
            call VecSetSizes(v%v, PETSC_DECIDE, vec_dat%nglob, ierr)
          else
            print *, "ERROR: invalid vector creation!"
            stop
          end if
        
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
  module subroutine axpy(alpha, x, y)

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
  module function norm(v, norm_type) result(n)

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

  !> @brief get pointer to data array 
  !
  !> param[in] vec - the vector object
  !> param[in/out] vec_data - the pointer to data array contained in vec
  !!
  !subroutine get_array_pointer(vec, vec_data)
  !  use petsc
  !  class(vector), intent(in) :: vec
  !  !real(accs_real), pointer, dimension(:) :: vec_data
  !  PetscScalar, pointer, dimension(:) :: vec_data

  !  select type (vec)
  !    type is (vector_petsc)
  !      print *, 'here 4'
  !      call VecGetArray(vec, vec_data)
  !      print *, 'here 5'
  !    
  !    class default
  !      print *, "Type unhandled"
  !      stop
  !  end select
  !end subroutine get_array_pointer

  module subroutine vec_view(vec)
    !use petscvec, only: VecView, PETSC_VIEWER_STDOUT_SELF
    !use petscvec, only: PetscViewer, PetscViewerBinaryOpen, PETSC_COMM_WORLD, FILE_MODE_WRITE
    use petscvec
    class(vector), intent(in) :: vec
    integer(accs_err) :: ierr
    type(tPetscViewer), pointer :: output_viewer

    !call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "solution.dat", FILE_MODE_WRITE, output_viewer)
    select type (vec)
      type is (vector_petsc)
        call VecView(vec%v, PETSC_VIEWER_STDOUT_SELF, ierr)
      class default
        print *, "Type unhandled 1"
        stop
    end select
  end subroutine vec_view

end submodule vec_petsc
