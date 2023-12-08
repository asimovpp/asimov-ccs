!v Program file for pressure-velocity coupling case
!
!  This case demonstrates solution of the Navier-Stokes equations
!  using the SIMPLE algorithm for pressure-velocity coupling.

program simple

  use petscvec
  use petscsys

  use ccs_base, only: mesh
  use constants, only: cell, face, field_u, field_v, field_w, field_p, field_p_prime, field_mf, &
                    field_viscosity, ccs_split_type_low_high, ccs_split_undefined
  use kinds, only: ccs_real, ccs_int
  use types, only: field, upwind_field, central_field, face_field, ccs_mesh, &
                   vector_spec, ccs_vector, fluid, fluid_solver_selector
  use parallel, only: initialise_parallel_environment, create_new_par_env, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments, sync
  use parallel_types, only: parallel_environment
  use mesh_utils, only: build_square_mesh
  use vec, only: create_vector, set_vector_location
  use petsctypes, only: vector_petsc
  use pv_coupling, only: solve_nonlinear
  use utils, only: set_size, initialise, update, get_field, add_field, &
                   get_fluid_solver_selector, set_fluid_solver_selector, &
                   allocate_fluid_fields
  use meshing, only: set_mesh_object, nullify_mesh_object

  implicit none

  class(parallel_environment), allocatable, target :: par_env
  class(parallel_environment), allocatable, target :: shared_env
  type(vector_spec) :: vec_sizes
  logical :: use_mpi_splitting

  class(field), allocatable :: u, v, w, p, p_prime, mf, viscosity, density

  integer(ccs_int) :: cps = 50 !< Default value for cells per side

  integer(ccs_int) :: it_start, it_end

  double precision :: start_time
  double precision :: end_time

  real(ccs_real) :: res_target = 1.0e-6 ! Default target residual

  logical :: u_sol = .true.  ! Solve u
  logical :: v_sol = .true.  ! Solve v
  logical :: w_sol = .false. ! Don't solve w
  logical :: p_sol = .true.  ! Solve p

  type(fluid) :: flow_fields
  type(fluid_solver_selector) :: fluid_sol

  ! Set start and end iteration numbers (eventually will be read from input file)
  it_start = 1
  it_end = 1000

  print *, "Starting SIMPLE demo"
  call initialise_parallel_environment(par_env)

  use_mpi_splitting = .false.
  call create_new_par_env(par_env, ccs_split_type_low_high, use_mpi_splitting, shared_env)

  call read_command_line_arguments(par_env)

  call sync(par_env)
  call timer(start_time)

  ! Create a square mesh
  print *, "Building mesh"
  mesh = build_square_mesh(par_env, shared_env, cps, 1.0_ccs_real)
  call set_mesh_object(mesh)

  ! Initialise fields
  print *, "Initialise fields"
  allocate (upwind_field :: u)
  allocate (upwind_field :: v)
  allocate (upwind_field :: w)
  allocate (central_field :: p)
  allocate (central_field :: p_prime)
  allocate (central_field :: viscosity)
  allocate (face_field :: mf)
  
  ! Create and initialise field vectors
  call initialise(vec_sizes)

  print *, "Create vectors"
  call set_vector_location(cell, vec_sizes)
  call set_size(par_env, mesh, vec_sizes)
  call create_vector(vec_sizes, u%values)
  call create_vector(vec_sizes, v%values)
  call create_vector(vec_sizes, w%values)
  call create_vector(vec_sizes, p%values)
  call create_vector(vec_sizes, p%x_gradients)
  call create_vector(vec_sizes, p%y_gradients)
  call create_vector(vec_sizes, p%z_gradients)
  call create_vector(vec_sizes, p_prime%values)
  call create_vector(vec_sizes, p_prime%x_gradients)
  call create_vector(vec_sizes, p_prime%y_gradients)
  call create_vector(vec_sizes, p_prime%z_gradients)
  call create_vector(vec_sizes, viscosity%values)
  call update(u%values)
  call update(v%values)
  call update(p%values)
  call update(p%x_gradients)
  call update(p%y_gradients)
  call update(p%z_gradients)
  call update(p_prime%values)
  call update(p_prime%x_gradients)
  call update(p_prime%y_gradients)
  call update(p_prime%z_gradients)
  call update(viscosity%values)

  call set_vector_location(face, vec_sizes)
  call set_size(par_env, mesh, vec_sizes)
  call create_vector(vec_sizes, mf%values)
  call update(mf%values)

  ! Initialise velocity field
  print *, "Initialise velocity field"
  call initialise_velocity(u, v, w, mf, viscosity)
  call update(u%values)
  call update(v%values)
  call update(mf%values)
  call update(viscosity%values)

  call set_fluid_solver_selector(field_u, u_sol, fluid_sol)
  call set_fluid_solver_selector(field_v, v_sol, fluid_sol)
  call set_fluid_solver_selector(field_w, w_sol, fluid_sol)
  call set_fluid_solver_selector(field_p, p_sol, fluid_sol)

  call add_field(u, flow_fields)
  call add_field(v, flow_fields)
  call add_field(w, flow_fields)
  call add_field(p, flow_fields)
  call add_field(p_prime, flow_fields)
  call add_field(mf, flow_fields)
  call add_field(viscosity, flow_fields) 
  call add_field(density, flow_fields) 

  ! Solve using SIMPLE algorithm
  print *, "Start SIMPLE"
  call solve_nonlinear(par_env, mesh, it_start, it_end, res_target, &
                       fluid_sol, flow_fields)

  ! Clean-up
  deallocate (u)
  deallocate (v)
  deallocate (w)
  deallocate (p)
  deallocate (p_prime)

  call timer(end_time)

  if (par_env%proc_id == par_env%root) then
    print *, "Elapsed time = ", (end_time - start_time)
  end if

  call nullify_mesh_object()
  call cleanup_parallel_environment(par_env)

contains

  subroutine initialise_velocity(u, v, w, mf, viscosity)

    use constants, only: add_mode
    use types, only: vector_values, cell_locator
    use meshing, only: create_cell_locator, get_global_index, get_local_num_cells
    use fv, only: calc_cell_coords
    use utils, only: set_values, set_mode, set_entry, set_row
    use vec, only: get_vector_data, restore_vector_data, create_vector_values

    ! Arguments
    class(field), intent(inout) :: u, v, w, mf, viscosity

    ! Local variables
    integer(ccs_int) :: n_local
    integer(ccs_int) :: row, col
    integer(ccs_int) :: local_idx, self_idx
    real(ccs_real) :: u_val, v_val, w_val
    type(cell_locator) :: self_loc
    type(vector_values) :: u_vals, v_vals, w_vals
    real(ccs_real), dimension(:), pointer :: u_data, v_data, w_data, mf_data, viscosity_data

    ! Set alias
    call get_local_num_cells(n_local)

    call create_vector_values(n_local, u_vals)
    call create_vector_values(n_local, v_vals)
    call create_vector_values(n_local, w_vals)

    call set_mode(add_mode, u_vals)
    call set_mode(add_mode, v_vals)
    call set_mode(add_mode, w_vals)

    ! Set initial values for velocity fields
    do local_idx = 1, n_local
      call create_cell_locator(local_idx, self_loc)
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

    call update(u%values)
    call update(v%values)
    call update(w%values)
    
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
    call get_vector_data(viscosity%values, viscosity_data)

    u_data(:) = 0.0_ccs_real
    v_data(:) = 0.0_ccs_real
    w_data(:) = 0.0_ccs_real
    mf_data(:) = 0.0_ccs_real
    viscosity_data(:) =  1.e-2_ccs_real

    call restore_vector_data(u%values, u_data)
    call restore_vector_data(v%values, v_data)
    call restore_vector_data(w%values, w_data)
    call restore_vector_data(mf%values, mf_data)
    call restore_vector_data(viscosity%values, viscosity_data)

  end subroutine initialise_velocity

end program simple
