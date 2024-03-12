!v Submodule file vec_petsc.smod
!
!  An implementation of vector objects using PETSc - the datatype and operations on it
!  (creation, destruction, setting/getting, ...)
!
!  @build petsc
submodule(vec) vec_petsc
#include "ccs_macros.inc"

  use kinds, only: ccs_err
  use petsctypes, only: vector_petsc
  use parallel_types_mpi, only: parallel_environment_mpi
  use constants, only: cell, face
  use petsc, only: ADD_VALUES, INSERT_VALUES, SCATTER_FORWARD
  use utils, only: debug_print, exit_print, str
  use error_codes

  implicit none

contains

  !> Create a PETSc-backed vector
  module subroutine create_vector(vec_properties, v, name)

    use petsc, only: PETSC_DECIDE, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE
    use petscvec, only: VecCreateGhost, VecSetSizes, VecSetFromOptions, VecSet, VecSetOption, &
                        VecCreate

    use meshing, only: get_local_num_cells, get_halo_num_cells, get_num_faces

    type(vector_spec), intent(in) :: vec_properties     !< the data describing how the vector should be created.
    class(ccs_vector), allocatable, intent(out) :: v    !< the vector specialised to type vector_petsc.
    character(len=*), optional, intent(in) :: name      !< name of the vector object

    integer(ccs_int), dimension(:), allocatable :: global_halo_indices
    integer(ccs_int) :: i
    integer(ccs_int) :: nlocal, nhalo
    integer(ccs_int) :: num_faces
    integer(ccs_err) :: ierr ! Error code

    allocate (vector_petsc :: v)

    select type (v)
    type is (vector_petsc)

      v%modeset = .false.
      v%checked_out = .false.
      if (present(name)) v%name = name

      select type (par_env => vec_properties%par_env)
      type is (parallel_environment_mpi)

        associate (mesh => vec_properties%mesh)

          select case (vec_properties%storage_location)
          case (cell)
            call get_local_num_cells(nlocal)
            call get_halo_num_cells(nhalo)
            associate (idx_global => mesh%topo%global_indices)
              allocate (global_halo_indices(nhalo))
              do i = 1, nhalo
                global_halo_indices(i) = idx_global(i + nlocal) - 1_ccs_int
              end do
              call VecCreateGhost(par_env%comm, &
                                  nlocal, PETSC_DECIDE, &
                                  nhalo, global_halo_indices, &
                                  v%v, ierr)
              deallocate (global_halo_indices)
            end associate
            ! Vector has ghost points, store this information
            v%ghosted = .true.
          case (face)
            call get_num_faces(num_faces)

            call VecCreate(par_env%comm, v%v, ierr)
            call VecSetSizes(v%v, num_faces, PETSC_DECIDE, ierr)

            ! Vector doesn't have ghost points, store this information
            v%ghosted = .false.
          end select
        end associate

        if (present(name)) then
          call VecSetOptionsPrefix(v%v, v%name // ':', ierr)
        end if
        call VecSetFromOptions(v%v, ierr)
        call VecSetOption(v%v, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE, ierr)
        call VecSet(v%v, 0.0_ccs_real, ierr)
        v%allocated = .true.

      class default
        call error_abort("Unknown parallel environment")

      end select

    class default
      call error_abort("Unknown vector type.")

    end select

  end subroutine

  !> Sets values in a PETSc vector
  module subroutine set_vector_values(val_dat, v)

    use petsc, only: VecSetValues

    use constants, only: insert_mode, add_mode

    class(*), intent(in) :: val_dat         !< contains the values, their indices and the mode to use for setting them.
    class(ccs_vector), intent(inout) :: v   !< the PETSc vector.

    integer(ccs_int) :: n    ! Number of elements to add
    integer(ccs_int) :: mode ! Append or insert mode
    integer(ccs_err) :: ierr ! Error code

    select type (v)
    type is (vector_petsc)

      select type (val_dat)
      type is (vector_values)

        ! First check if safe to set
        if (v%modeset) then
          if (val_dat%setter_mode /= v%mode) then
            print *, ""
            call error_abort("ERROR: trying to set vector using different mode without updating.")
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
          call error_abort("Unknown mode.")
        end if

        call VecSetValues(v%v, n, val_dat%global_indices, val_dat%values, mode, ierr)

      class default
        call error_abort("Unknown vector value type.")

      end select

    class default
      call error_abort("Unknown vector type.")

    end select

  end subroutine

  !> Perform a parallel update of a PETSc vector
  module subroutine update_vector(v)

    class(ccs_vector), intent(inout) :: v   !< the PETSc vector

    select type (v)
    type is (vector_petsc)

      call begin_update_vector(v)
      call end_update_vector(v)

      call begin_ghost_update_vector(v)
      call end_ghost_update_vector(v)

    class default
      call error_abort("Unknown vector type.")

    end select

  end subroutine

  !v Begin a parallel update of a PETSc vector
  !
  !  Begins the parallel update to allow overlapping comms and compute
  module subroutine begin_update_vector(v)

    use petsc, only: VecAssemblyBegin

    class(ccs_vector), intent(inout) :: v   !< the PETSc vector

    integer(ccs_err) :: ierr ! Error code

    select type (v)
    type is (vector_petsc)

      call VecAssemblyBegin(v%v, ierr)

    class default
      call error_abort("Unknown vector type.")

    end select

  end subroutine

  !v End a parallel update of a PETSc vector.
  !
  !  Ends the parallel update to allow overlapping comms and compute.
  module subroutine end_update_vector(v)

    use petsc, only: VecAssemblyEnd

    class(ccs_vector), intent(inout) :: v   !< the PETSc vector

    integer(ccs_err) :: ierr ! Error code

    select type (v)
    type is (vector_petsc)

      call VecAssemblyEnd(v%v, ierr)

      v%modeset = .false. ! It is now safe to change value setting mode

    class default
      call error_abort("Unknown vector type.")

    end select

  end subroutine

  !v Begin a ghost update of a PETSc vector
  !
  !  Begins the ghost update to allow overlapping comms and compute
  subroutine begin_ghost_update_vector(v)

    use petsc, only: VecGhostUpdateBegin

    class(ccs_vector), intent(inout) :: v   !< the PETSc vector

    integer(ccs_err) :: ierr ! Error code

    select type (v)
    type is (vector_petsc)

      if (v%ghosted) then
        ! Cant update ghosts if not ghost points.
        call VecGhostUpdateBegin(v%v, INSERT_VALUES, SCATTER_FORWARD, ierr)
      end if

    class default
      call error_abort("Unknown vector type.")

    end select

  end subroutine

  !v End a ghost update of a PETSc vector.
  !
  !  Ends the ghost update to allow overlapping comms and compute.
  subroutine end_ghost_update_vector(v)

    use petsc, only: VecGhostUpdateEnd

    class(ccs_vector), intent(inout) :: v   !< the PETSc vector

    integer(ccs_err) :: ierr ! Error code

    select type (v)
    type is (vector_petsc)

      if (v%ghosted) then
        ! Cant update ghosts if not ghost points.
        call VecGhostUpdateEnd(v%v, INSERT_VALUES, SCATTER_FORWARD, ierr)
      end if

    class default
      call error_abort("Unknown vector type.")

    end select

  end subroutine

  !v Perform the AXPY vector operation using PETSc
  !
  !          y[i] = alpha * x[i] + y[i]
  module subroutine vec_axpy(alpha, x, y)

    use petscvec, only: VecAXPY

    real(ccs_real), intent(in) :: alpha     !< a scalar value
    class(ccs_vector), intent(in) :: x      !< a PETSc input vector
    class(ccs_vector), intent(inout) :: y   !< PETSc vector serving as input, overwritten with result

    integer(ccs_err) :: ierr ! Error code

    select type (x)
    type is (vector_petsc)

      select type (y)
      type is (vector_petsc)

        ! PETSc performs AXPY as YPAX, with result stored in Y.
        call VecAXPY(y%v, alpha, x%v, ierr)

      class default
        call error_abort("Unknown vector type.")

      end select

    class default
      call error_abort("Unknown vector type.")

    end select

  end subroutine

  !v Perform the AYPX vector operation using PETSc
  !
  !          y[i] = x[i] + beta * y[i]
  module subroutine vec_aypx(x, beta, y)

    use petscvec, only: VecAYPX

    real(ccs_real), intent(in) :: beta      !< a scalar value
    class(ccs_vector), intent(in) :: x      !< a PETSc input vector
    class(ccs_vector), intent(inout) :: y   !< PETSc vector serving as input, overwritten with result

    integer(ccs_err) :: ierr ! Error code

    select type (x)
    type is (vector_petsc)

      select type (y)
      type is (vector_petsc)

        call VecAYPX(y%v, beta, x%v, ierr)

      class default
        call error_abort("Unknown vector type.")

      end select

    class default
      call error_abort("Unknown vector type.")

    end select

  end subroutine

  !> Compute the norm of a PETSc vector
  module function vec_norm(v, norm_type) result(n)

    use petscvec, only: NORM_2, NORM_INFINITY, VecNorm

    class(ccs_vector), intent(in) :: v          !< the PETSc vector
    integer(ccs_int), intent(in) :: norm_type   !< which norm to compute? Currently supported is the 2 norm: norm_type=2.

    real(ccs_real) :: n      !< The computed norm
    integer(ccs_err) :: ierr ! Error code

    n = 0.0_ccs_real ! initialise norm to 0

    select type (v)
    type is (vector_petsc)

      select case (norm_type)
      case (0)
        call VecNorm(v%v, NORM_INFINITY, n, ierr)
      case (2)
        call VecNorm(v%v, NORM_2, n, ierr)
      case default
        call error_abort("ERROR: unknown vector norm type " // str(norm_type))
      end select

    class default
      call error_abort("Type unhandled")
    end select

  end function

  pure module subroutine clear_vector_values_entries(val_dat)
    type(vector_values), intent(inout) :: val_dat

    val_dat%global_indices(:) = -1 ! PETSc ignores -ve indices, used as "empty" indicator
    val_dat%values(:) = 0.0_ccs_real

  end subroutine clear_vector_values_entries

  pure module subroutine set_vector_values_row(row, val_dat)
    integer(ccs_int), intent(in) :: row
    type(vector_values), intent(inout) :: val_dat

    integer(ccs_int), dimension(1) :: idxs !< Temporary array mapping rows to indices in the
    !< current working set. N.B. the dimension of this
    !< array must match the rank of
    !< vector_values%global_indices.
    integer(ccs_int) :: i
    integer(ccs_int) :: petsc_row

    petsc_row = row - 1 ! PETSc is zero-indexed

    idxs = findloc(val_dat%global_indices, petsc_row, kind=ccs_int)
    i = idxs(1) ! We want the first entry
    if (i == 0) then
      ! New entry
      idxs = findloc(val_dat%global_indices, -1_ccs_int, kind=ccs_int)
      i = idxs(1) ! We want the first entry
      if (i == 0) then
        error stop no_free_entry ! Couldn't find a free entry in vector values
      end if
    end if

    val_dat%current_entry = i
    val_dat%global_indices(i) = petsc_row

  end subroutine set_vector_values_row

  !> Gets the data in a given vector with possibility to overwrite the data.
  module subroutine get_vector_data(vec, array)
    use petscvec, only: VecGhostGetLocalForm, VecGetArrayF90
    class(ccs_vector), intent(inout) :: vec !< the vector to get data from
    real(ccs_real), dimension(:), pointer, intent(out) :: array !< an array to store the data in
    integer :: ierr

    select type (vec)
    type is (vector_petsc)
      if (vec%checked_out) then
        call error_abort("ERROR: trying to access already checked-out vector")
      end if

      if (vec%modeset) then
        call error_abort("WARNING: trying to access vector without updating")
      end if

      if (vec%ghosted) then
        call VecGhostGetLocalForm(vec%v, vec%v_local, ierr)
        call VecGetArrayF90(vec%v_local, array, ierr)
      else
        call VecGetArrayF90(vec%v, array, ierr)
      end if

      vec%checked_out = .true.
    class default
      call error_abort('Invalid vector type.')
    end select
  end subroutine get_vector_data

  !> Resets the vector data if required for further processing
  module subroutine restore_vector_data(vec, array)
    use petscvec, only: VecRestoreArrayF90, VecGhostRestoreLocalForm

    class(ccs_vector), intent(inout) :: vec !< the vector to reset
    real(ccs_real), dimension(:), pointer, intent(in) :: array !< the array containing the data to restore

    integer :: ierr

    select type (vec)
    type is (vector_petsc)
      if (.not. vec%checked_out) then
        call error_abort("ERROR: trying to double-restore vector")
      end if

      if (vec%ghosted) then
        call VecRestoreArrayF90(vec%v_local, array, ierr)
        call VecGhostRestoreLocalForm(vec%v, vec%v_local, ierr)
      else
        call VecRestoreArrayF90(vec%v, array, ierr)
      end if

      vec%checked_out = .false.

    class default
      call error_abort('Invalid vector type.')
    end select
  end subroutine restore_vector_data

  !> Gets the data in a given vector with readonly access.
  module subroutine get_vector_data_readonly(vec, array)
    use petscvec, only: VecGhostGetLocalForm, VecGetArrayReadF90
    class(ccs_vector), intent(inout) :: vec !< the vector to get data from
    real(ccs_real), dimension(:), pointer, intent(out) :: array !< an array to store the data in
    integer :: ierr

    select type (vec)
    type is (vector_petsc)
      if (vec%ghosted) then
        call VecGhostGetLocalForm(vec%v, vec%v_local, ierr)
        call VecGetArrayReadF90(vec%v_local, array, ierr)
      else
        call VecGetArrayReadF90(vec%v, array, ierr)
      end if
    class default
      call error_abort('Invalid vector type.')
    end select
  end subroutine get_vector_data_readonly

  !> Resets the vector data with readonly access.
  module subroutine restore_vector_data_readonly(vec, array)
    use petscvec, only: VecRestoreArrayReadF90, VecGhostRestoreLocalForm

    class(ccs_vector), intent(inout) :: vec !< the vector to reset
    real(ccs_real), dimension(:), pointer, intent(in) :: array !< the array containing the data to restore

    integer :: ierr

    select type (vec)
    type is (vector_petsc)
      if (vec%ghosted) then
        call VecRestoreArrayReadF90(vec%v_local, array, ierr)
        call VecGhostRestoreLocalForm(vec%v, vec%v_local, ierr)
      else
        call VecRestoreArrayReadF90(vec%v, array, ierr)
      end if
    class default
      call error_abort('Invalid vector type.')
    end select
  end subroutine restore_vector_data_readonly

  module subroutine zero_vector(vec)
    use petscvec, only: VecZeroEntries

    class(ccs_vector), intent(inout) :: vec

    integer(ccs_err) :: ierr

    select type (vec)
    type is (vector_petsc)
      call VecZeroEntries(vec%v, ierr)
    class default
      call error_abort('Invalid vector type.')
    end select

  end subroutine zero_vector

  !> Replaces each component of a vector by its reciprocal
  module subroutine vec_reciprocal(vec)

    use petscvec, only: VecReciprocal

    class(ccs_vector), intent(inout) :: vec !< the vector

    integer(ccs_err) :: ierr ! Error code

    select type (vec)
    type is (vector_petsc)

      call VecReciprocal(vec%v, ierr)

    class default
      call error_abort("Unknown vector type.")
    end select
  end subroutine vec_reciprocal

  module subroutine mult_vec_vec(a, b)
    use petscvec, only: VecPointwiseMult

    class(ccs_vector), intent(in) :: a
    class(ccs_vector), intent(inout) :: b

    integer(ccs_err) :: ierr

    select type (a)
    type is (vector_petsc)
      select type (b)
      type is (vector_petsc)
        call VecPointwiseMult(b%v, a%v, b%v, ierr)
      class default
        call error_abort("Unknown vector type.")
      end select
    class default
      call error_abort("Unknown vector type.")
    end select

  end subroutine mult_vec_vec

  module subroutine scale_vec(alpha, v)
    use petscvec, only: VecScale

    real(ccs_real), intent(in) :: alpha
    class(ccs_vector), intent(inout) :: v

    integer(ccs_err) :: ierr

    select type (v)
    type is (vector_petsc)
      call VecScale(v%v, alpha, ierr)
    class default
      call error_abort("Unknown vector type.")
    end select

  end subroutine

  !> PETSc implementation to get vector data in natural ordering
  module subroutine get_natural_data_vec(par_env, mesh, v, data)

    use petsc, only: PETSC_DECIDE
    use petscvec, only: tVec, &
                        VecCreate, VecSetSizes, VecSetFromOptions, &
                        VecSetValues, VecAssemblyBegin, VecAssemblyEnd, &
                        VecGetArrayReadF90, VecRestoreArrayReadF90, &
                        VecDestroy

    class(parallel_environment), intent(in) :: par_env
    type(ccs_mesh), intent(in) :: mesh
    class(ccs_vector), intent(inout) :: v
    real(ccs_real), dimension(:), allocatable, intent(out) :: data !< The returned vector data in
    !< natural ordering. Note the use
    !< of allocatable + intent(out),
    !< this ensures it will be
    !< de/reallocated by this subroutine.

    real(ccs_real), dimension(:), pointer :: vec_data ! The data stored in the vector

    type(tVec) :: vec_tmp

    integer(ccs_err) :: ierr

    associate (topo => mesh%topo, &
               local_num_cells => mesh%topo%local_num_cells)

      ! Create temporary vector to hold the shuffled data
      select type (par_env)
      type is (parallel_environment_mpi)
        call VecCreate(par_env%comm, vec_tmp, ierr)
      class default
        call error_abort("Only MPI parallel environments currently supported!")
      end select
      call VecSetSizes(vec_tmp, local_num_cells, PETSC_DECIDE, ierr)
      call VecSetFromOptions(vec_tmp, ierr)

      ! Copy local data into temporary vector, inserting values at natural index locations
      call get_vector_data(v, vec_data)
      call VecSetValues(vec_tmp, local_num_cells, topo%natural_indices(1:local_num_cells) - 1, &
                        vec_data(1:local_num_cells), INSERT_VALUES, ierr)
      call VecAssemblyBegin(vec_tmp, ierr)
      call VecAssemblyEnd(vec_tmp, ierr)
      call restore_vector_data(v, vec_data)

      ! Extract natural-indexed data into data array and return
      call VecGetArrayReadF90(vec_tmp, vec_data, ierr)
      if (allocated(data)) then ! Really shouldn't happen
        deallocate (data)
      end if
      allocate (data(local_num_cells))
      data(1:local_num_cells) = vec_data(1:local_num_cells)
      call VecRestoreArrayReadF90(vec_tmp, vec_data, ierr)

      ! Clean up
      call VecDestroy(vec_tmp, ierr)

    end associate

  end subroutine get_natural_data_vec

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
!           call error_abort("Unknown parallel environment")
!       end select

!       class default
!         call error_abort("Unknown vector type.")
!     end select
!   end subroutine vec_view

end submodule vec_petsc
