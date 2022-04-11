!> @brief Program file for pressure-velocity coupling case
!
!> @details This case demonstrates solution of the Navier-Stokes equations
!!          using the SIMPLE algorithm for pressure-velocity coupling.
!

program simple

  use petscvec
  use petscsys

  use constants, only : cell, face
  use kinds, only: accs_real, accs_int
  use types, only: field, upwind_field, central_field, face_field, mesh, &
                   vector_init_data, matrix, vector
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments, sync
  use parallel_types, only: parallel_environment
  use mesh_utils, only: build_square_mesh
  use vec, only: create_vector, set_vector_location
  use petsctypes, only: matrix_petsc, vector_petsc
  use pv_coupling, only: solve_nonlinear
  use utils, only: set_size, initialise, update
                      
  implicit none

  class(parallel_environment), allocatable, target :: par_env
  type(mesh)             :: square_mesh
  type(vector_init_data) :: vec_sizes

  class(field), allocatable :: u, v, p, pp, mf

  integer(accs_int) :: cps = 50 ! Default value for cells per side

  integer(accs_int) :: it_start, it_end, ierr

  double precision :: start_time
  double precision :: end_time
 
  type(tPetscViewer) :: viewer

  ! Set start and end iteration numbers (eventually will be read from input file)
  it_start = 1
  it_end   = 1000

  print *, "Starting SIMPLE demo"
  call initialise_parallel_environment(par_env)
  call read_command_line_arguments(par_env)

  call sync(par_env)
  call timer(start_time)

  ! Create a square mesh
  print *, "Building mesh"
  square_mesh = build_square_mesh(par_env, cps, 1.0_accs_real)

  ! Initialise fields
  print *, "Initialise fields"
  allocate(upwind_field :: u)
  allocate(upwind_field :: v)
  allocate(central_field :: p)
  allocate(central_field :: pp)
  allocate(face_field :: mf)

  ! Create and initialise field vectors
  call initialise(vec_sizes)

  print *, "Create vectors"
  call set_vector_location(cell, vec_sizes)
  call set_size(par_env, square_mesh, vec_sizes)
  call create_vector(vec_sizes, u%vec)
  call create_vector(vec_sizes, v%vec)
  call create_vector(vec_sizes, p%vec)
  call create_vector(vec_sizes, p%gradx)
  call create_vector(vec_sizes, p%grady)
  call create_vector(vec_sizes, p%gradz)
  call create_vector(vec_sizes, pp%vec)
  call create_vector(vec_sizes, pp%gradx)
  call create_vector(vec_sizes, pp%grady)
  call create_vector(vec_sizes, pp%gradz)
  call update(u%vec)
  call update(v%vec)
  call update(p%vec)
  call update(p%gradx)
  call update(p%grady)
  call update(p%gradz)
  call update(pp%vec)
  call update(pp%gradx)
  call update(pp%grady)
  call update(pp%gradz)

  call set_vector_location(face, vec_sizes)
  call set_size(par_env, square_mesh, vec_sizes)
  call create_vector(vec_sizes, mf%vec)
  call update(mf%vec)
  
  ! Initialise velocity field
  print *, "Initialise velocity field"
  call initialise_velocity(square_mesh, u, v, mf)
  call update(u%vec)
  call update(v%vec)
  call update(mf%vec)

  ! Solve using SIMPLE algorithm
  print *, "Start SIMPLE"
  call solve_nonlinear(par_env, square_mesh, cps, it_start, it_end, u, v, p, pp, mf)

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,"u",FILE_MODE_WRITE,viewer, ierr)

  associate (vec => u%vec)
    select type (vec)
    type is (vector_petsc)
      call VecView(vec%v, viewer, ierr)
    end select
  end associate

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,"v",FILE_MODE_WRITE,viewer, ierr)

  associate (vec => v%vec)
    select type (vec)
    type is (vector_petsc)
      call VecView(vec%v, viewer, ierr)
    end select
  end associate

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,"p",FILE_MODE_WRITE,viewer, ierr)

  associate (vec => p%vec)
    select type (vec)
    type is (vector_petsc)
      call VecView(vec%v, viewer, ierr)
    end select
  end associate

  call PetscViewerDestroy(viewer,ierr)

  ! Clean-up
  deallocate(u)
  deallocate(v)
  deallocate(p)
  deallocate(pp)

  call timer(end_time)

  if (par_env%proc_id == par_env%root) then
    print *, "Elapsed time = ", (end_time - start_time)
  end if

  call cleanup_parallel_environment(par_env)

contains

  subroutine initialise_velocity(cell_mesh, u, v, mf)

    use constants, only: add_mode
    use types, only: vector_values, cell_locator
    use meshing, only: set_cell_location, get_global_index
    use fv, only: calc_cell_coords
    use utils, only: pack_entries, set_values
    use vec, only : get_vector_data, restore_vector_data
    
    ! Arguments
    class(mesh), intent(in) :: cell_mesh
    class(field), intent(inout) :: u, v, mf

    ! Local variables
    integer(accs_int) :: row, col
    integer(accs_int) :: local_idx, self_idx
    real(accs_real) :: u_val, v_val
    type(cell_locator) :: self_loc
    type(vector_values) :: u_vals, v_vals
    real(accs_real), dimension(:), pointer :: u_data, v_data, mf_data

    ! Set mode
    u_vals%mode = add_mode
    v_vals%mode = add_mode

    ! Set alias
    associate(n_local => cell_mesh%nlocal)
      ! Allocate temporary arrays for storing global cell indices 
      allocate(u_vals%idx(n_local))
      allocate(v_vals%idx(n_local))

      ! Allocate temporary arrays for storing values
      allocate(u_vals%val(n_local))
      allocate(v_vals%val(n_local))

      ! Set initial values for velocity fields
      do local_idx = 1, n_local
        call set_cell_location(cell_mesh, local_idx, self_loc)
        call get_global_index(self_loc, self_idx)
        call calc_cell_coords(self_idx, cps, row, col)

        u_val = real(col, accs_real)/real(cps, accs_real)
        v_val = -real(row, accs_real)/real(cps, accs_real)

        call pack_entries(local_idx, self_idx, u_val, u_vals)
        call pack_entries(local_idx, self_idx, v_val, v_vals)
      end do
    end associate

    call set_values(u_vals, u%vec)
    call set_values(v_vals, v%vec)

    deallocate(u_vals%idx, v_vals%idx, u_vals%val, v_vals%val)

    call get_vector_data(u%vec, u_data)
    call get_vector_data(v%vec, v_data)
    call get_vector_data(mf%vec, mf_data)

    u_data(:) = 0.0_accs_real
    v_data(:) = 0.0_accs_real
    mf_data(:) = 0.0_accs_real
    
    call restore_vector_data(u%vec, u_data)
    call restore_vector_data(v%vec, v_data)
    call restore_vector_data(mf%vec, mf_data)
    
  end subroutine initialise_velocity

end program simple
