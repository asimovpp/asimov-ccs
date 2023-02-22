!v Program file for pressure-velocity coupling case
!
!  This case demonstrates solution of the Navier-Stokes equations
!  using the SIMPLE algorithm for pressure-velocity coupling.

program simple

  use petscvec
  use petscsys

  use constants, only: cell, face
  use kinds, only: ccs_real, ccs_int
  use types, only: field, upwind_field, central_field, face_field, ccs_mesh, &
                   vector_spec, ccs_vector
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments, sync
  use parallel_types, only: parallel_environment
  use mesh_utils, only: build_square_mesh
  use vec, only: create_vector, set_vector_location
  use petsctypes, only: vector_petsc
  use pv_coupling, only: solve_nonlinear
  use utils, only: set_size, initialise, update

  implicit none

  class(parallel_environment), allocatable, target :: par_env
  type(ccs_mesh) :: square_mesh
  type(vector_spec) :: vec_sizes

  class(field), allocatable :: u, v, w, p, pp, mf

  integer(ccs_int) :: cps = 50 !< Default value for cells per side

  integer(ccs_int) :: it_start, it_end

  double precision :: start_time
  double precision :: end_time

  real(ccs_real) :: res_target = 1.0e-6 ! Default target residual

  logical :: u_sol = .true.  ! Solve u
  logical :: v_sol = .true.  ! Solve v
  logical :: w_sol = .false. ! Don't solve w
  logical :: p_sol = .true.  ! Solve p

  ! Set start and end iteration numbers (eventually will be read from input file)
  it_start = 1
  it_end = 1000

  print *, "Starting SIMPLE demo"
  call initialise_parallel_environment(par_env)
  call read_command_line_arguments(par_env)

  call sync(par_env)
  call timer(start_time)

  ! Create a square mesh
  print *, "Building mesh"
  square_mesh = build_square_mesh(par_env, cps, 1.0_ccs_real)

  ! Initialise fields
  print *, "Initialise fields"
  allocate (upwind_field :: u)
  allocate (upwind_field :: v)
  allocate (upwind_field :: w)
  allocate (central_field :: p)
  allocate (central_field :: pp)
  allocate (face_field :: mf)

  ! Create and initialise field vectors
  call initialise(vec_sizes)

  print *, "Create vectors"
  call set_vector_location(cell, vec_sizes)
  call set_size(par_env, square_mesh, vec_sizes)
  call create_vector(vec_sizes, u%values)
  call create_vector(vec_sizes, v%values)
  call create_vector(vec_sizes, w%values)
  call create_vector(vec_sizes, p%values)
  call create_vector(vec_sizes, p%x_gradients)
  call create_vector(vec_sizes, p%y_gradients)
  call create_vector(vec_sizes, p%z_gradients)
  call create_vector(vec_sizes, pp%values)
  call create_vector(vec_sizes, pp%x_gradients)
  call create_vector(vec_sizes, pp%y_gradients)
  call create_vector(vec_sizes, pp%z_gradients)
  call update(u%values)
  call update(v%values)
  call update(p%values)
  call update(p%x_gradients)
  call update(p%y_gradients)
  call update(p%z_gradients)
  call update(pp%values)
  call update(pp%x_gradients)
  call update(pp%y_gradients)
  call update(pp%z_gradients)

  call set_vector_location(face, vec_sizes)
  call set_size(par_env, square_mesh, vec_sizes)
  call create_vector(vec_sizes, mf%values)
  call update(mf%values)

  ! Initialise velocity field
  print *, "Initialise velocity field"
  call initialise_velocity(square_mesh, u, v, w, mf)
  call update(u%values)
  call update(v%values)
  call update(mf%values)

  ! Solve using SIMPLE algorithm
  print *, "Start SIMPLE"
  call solve_nonlinear(par_env, square_mesh, it_start, it_end, res_target, &
                       u_sol, v_sol, w_sol, p_sol, u, v, w, p, pp, mf)

  ! Clean-up
  deallocate (u)
  deallocate (v)
  deallocate (w)
  deallocate (p)
  deallocate (pp)

  call timer(end_time)

  if (par_env%proc_id == par_env%root) then
    print *, "Elapsed time = ", (end_time - start_time)
  end if

  call cleanup_parallel_environment(par_env)

contains

  subroutine initialise_velocity(cell_mesh, u, v, w, mf)

    use constants, only: add_mode
    use types, only: vector_values, cell_locator
    use meshing, only: set_cell_location, get_global_index, get_local_num_cells
    use fv, only: calc_cell_coords
    use utils, only: set_values, set_mode, set_entry, set_row
    use vec, only: get_vector_data, restore_vector_data, create_vector_values

    ! Arguments
    class(ccs_mesh), intent(in) :: cell_mesh
    class(field), intent(inout) :: u, v, w, mf

    ! Local variables
    integer(ccs_int) :: n_local
    integer(ccs_int) :: row, col
    integer(ccs_int) :: local_idx, self_idx
    real(ccs_real) :: u_val, v_val, w_val
    type(cell_locator) :: self_loc
    type(vector_values) :: u_vals, v_vals, w_vals
    real(ccs_real), dimension(:), pointer :: u_data, v_data, w_data, mf_data

    ! Set alias
    call get_local_num_cells(cell_mesh, n_local)

    call create_vector_values(n_local, u_vals)
    call create_vector_values(n_local, v_vals)
    call create_vector_values(n_local, w_vals)

    call set_mode(add_mode, u_vals)
    call set_mode(add_mode, v_vals)
    call set_mode(add_mode, w_vals)

    ! Set initial values for velocity fields
    do local_idx = 1, n_local
       call set_cell_location(cell_mesh, local_idx, self_loc)
       call get_global_index(self_loc, self_idx)
       call calc_cell_coords(self_idx, cps, row, col)

       u_val = real(col, ccs_real) / real(cps, ccs_real)
       v_val = -real(row, ccs_real) / real(cps, ccs_real)
       w_val = 0.0_ccs_real

       call set_row(self_idx, u_vals)
       call set_entry(u_val, u_vals)

       call set_row(self_idx, v_vals)
       call set_entry(v_val, v_vals)

       call set_row(self_idx, w_vals)
       call set_entry(w_val, w_vals)
    end do

    call set_values(u_vals, u%values)
    call set_values(v_vals, v%values)
    call set_values(w_vals, w%values)

    ! XXX: make a finaliser for vector values
    deallocate (u_vals%global_indices)
    deallocate (v_vals%global_indices)
    deallocate (w_vals%global_indices)
    deallocate (u_vals%values)
    deallocate (v_vals%values)
    deallocate (w_vals%values)

    call get_vector_data(u%values, u_data)
    call get_vector_data(v%values, v_data)
    call get_vector_data(w%values, w_data)
    call get_vector_data(mf%values, mf_data)

    u_data(:) = 0.0_ccs_real
    v_data(:) = 0.0_ccs_real
    w_data(:) = 0.0_ccs_real
    mf_data(:) = 0.0_ccs_real

    call restore_vector_data(u%values, u_data)
    call restore_vector_data(v%values, v_data)
    call restore_vector_data(w%values, w_data)
    call restore_vector_data(mf%values, mf_data)

  end subroutine initialise_velocity

end program simple
