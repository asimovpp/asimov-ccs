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


  module procedure clear_vector_values_entries

    val_dat%idx(:) = -1 ! PETSc ignores -ve indices, used as "empty" indicator
    val_dat%val(:) = 0.0_accs_real
    
  end procedure clear_vector_values_entries
  
  module procedure set_vector_values_row

    integer(accs_int), dimension(rank(val_dat%idx)) :: idxs
    integer(accs_int) :: i
    logical :: new_entry
    integer(accs_int) :: petsc_row

    petsc_row = row - 1 ! PETSc is zero-indexed
    new_entry = .false.
    
    idxs = findloc(val_dat%idx, petsc_row, kind=accs_int)
    i = idxs(1) ! We want the first entry
    if (i == 0) then
      new_entry = .true.
    end if

    if (new_entry) then
      idxs = findloc(val_dat%idx, -1_accs_int, kind=accs_int)
      i = idxs(1) ! We want the first entry
      if (i == 0) then
        print *, "ERROR: Couldn't find a free entry in vector values!"
        stop
      end if
    end if
    
    val_dat%current_entry = i
    val_dat%idx(i) = petsc_row
    
  end procedure set_vector_values_row
  
end submodule vec_petsc
