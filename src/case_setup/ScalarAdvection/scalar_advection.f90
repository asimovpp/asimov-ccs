!v Program file for scalar advection case

program scalar_advection
#include "ccs_macros.inc"

  ! ASiMoV-CCS uses
  use ccs_base, only: mesh
  use kinds, only: ccs_real, ccs_int
  use case_config, only: num_steps, num_iters, dt, cps, domain_size, write_frequency, &
                         velocity_relax, pressure_relax, res_target, case_name, write_gradients, &
                         velocity_solver_method_name, velocity_solver_precon_name, &
                         pressure_solver_method_name, pressure_solver_precon_name, restart, unsteady
  use types, only: vector_spec, ccs_vector, matrix_spec, ccs_matrix, field_spec, &
                   equation_system, linear_solver, ccs_mesh, field_ptr, &
                   field, upwind_field, central_field, bc_config, face_locator, &
                   fluid
  use constants, only: cell, face, ccsconfig, ccs_string_len, geoext, adiosconfig, ndim, &
                       cell_centred_central, cell_centred_upwind, face_centred, &
                       ccs_split_type_shared, ccs_split_type_low_high, ccs_split_undefined
  use meshing, only: get_boundary_status, create_face_locator, get_total_num_cells, get_global_num_cells
  use fields, only: create_field, set_field_config_file, set_field_n_boundaries, set_field_name, &
       set_field_type, set_field_vector_properties, set_field_enable_cell_corrections
  use fortran_yaml_c_interface, only: parse
  use vec, only: create_vector, set_vector_location
  use mat, only: create_matrix, set_nnz
  use solver, only: create_solver, solve, set_equation_system
  use utils, only: update, initialise, set_size, add_field_to_outputlist, exit_print, finalise, zero, &
                   get_field
  use mesh_utils, only: build_square_mesh, write_mesh, compute_face_interpolation
  use meshing, only: set_mesh_object, nullify_mesh_object
  use parallel_types, only: parallel_environment
  use parallel, only: initialise_parallel_environment, create_new_par_env, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments, is_root
  use fv, only: compute_fluxes, update_gradient
  use io_visualisation, only: write_solution, read_solution
  use read_config, only: get_variables, get_boundary_count, get_boundary_names, get_case_name, &
                         get_enable_cell_corrections, get_variable_types
  use boundary_conditions, only: read_bc_config, allocate_bc_arrays
  use timers, only: timer_init, timer_register_start, timer_register, timer_start, timer_stop, timer_print, timer_get_time, timer_print_all

  implicit none

  class(parallel_environment), allocatable, target :: par_env
  class(parallel_environment), allocatable, target :: shared_env
  class(ccs_vector), allocatable, target :: source
  class(ccs_matrix), allocatable, target :: M
  class(linear_solver), allocatable :: scalar_solver
  character(len=:), allocatable :: input_path  ! Path to input directory
  character(len=:), allocatable :: case_path  ! Path to input directory with case name appended
  character(len=:), allocatable :: ccs_config_file ! Config file for CCS
  character(len=ccs_string_len), dimension(:), allocatable :: variable_names  ! variable names for BC reading
  integer(ccs_int), dimension(:), allocatable :: variable_types              ! cell centred upwind, central, etc.

  type(vector_spec) :: vec_properties
  type(matrix_spec) :: mat_properties
  type(equation_system) :: scalar_equation_system

  integer(ccs_int) :: n_boundaries
  logical :: enable_cell_corrections

  type(field_spec) :: field_properties
  class(field), pointer :: u, v, mf, viscosity, density
  class(field), pointer :: scalar

  integer(ccs_int) :: direction = 0 ! pass zero for "direction" of scalar field when computing fluxes
  integer(ccs_int) :: irank ! MPI rank ID
  integer(ccs_int) :: isize ! Size of MPI world
  integer(ccs_int) :: i

  double precision :: start_time
  double precision :: end_time
  logical :: use_mpi_splitting

  type(fluid) :: flow_fields

  character(len=128), dimension(:), allocatable :: bnd_names
  integer(ccs_int):: t          ! Timestep counter
  integer(ccs_int):: timer_index_sol
  integer(ccs_int):: timer_index_io_sol

  call initialise_parallel_environment(par_env)
  use_mpi_splitting = .false.
  call create_new_par_env(par_env, ccs_split_type_low_high, use_mpi_splitting, shared_env)

  call timer_init()

  irank = par_env%proc_id
  isize = par_env%num_procs

  call read_command_line_arguments(par_env, cps, case_name=case_name, in_dir=input_path)

  if (allocated(input_path)) then
    case_path = input_path // "/" // case_name
  else
    case_path = case_name
  end if

  ccs_config_file = case_path // ccsconfig

  call timer(start_time)

  ! Read case name and runtime parameters from configuration file
  call read_configuration(ccs_config_file)

  if (irank == par_env%root) print *, "Starting ", case_name, " case!"

  ! set solver and preconditioner info
  velocity_solver_method_name = "gmres"
  velocity_solver_precon_name = "bjacobi"
  pressure_solver_method_name = "cg"
  pressure_solver_precon_name = "gamg"

  ! Set up the square mesh
  if (irank == par_env%root) print *, "Building mesh"
  call get_boundary_names(ccs_config_file, bnd_names)
  mesh = build_square_mesh(par_env, shared_env, cps, 1.0_ccs_real, bnd_names)
  call set_mesh_object(mesh)

  ! Initialise fields
  if (irank == par_env%root) print *, "Initialise fields"

  ! Write gradients to solution file
  write_gradients = .true.

  ! Read boundary conditions
  if (irank == par_env%root) print *, "Read and allocate BCs"
  call get_boundary_count(ccs_config_file, n_boundaries)

  ! Create and initialise field vectors
  if (irank == par_env%root) print *, "Initialise field vectors"
  call initialise(vec_properties)

  call set_vector_location(cell, vec_properties)
  call set_size(par_env, mesh, vec_properties)
  call get_enable_cell_corrections(ccs_config_file, enable_cell_corrections)

  call set_field_config_file(ccs_config_file, field_properties)
  call set_field_n_boundaries(n_boundaries, field_properties)
  call set_field_enable_cell_corrections(enable_cell_corrections, field_properties)

  call set_field_vector_properties(vec_properties, field_properties)

  if (is_root(par_env)) then
    print *, "Build field list"
  end if

  do i = 1, size(variable_names)
    if (is_root(par_env)) then
      print *, "Creating field ", trim(variable_names(i))
    end if
    call set_field_type(variable_types(i), field_properties)
    call set_field_name(variable_names(i), field_properties)
    call create_field(par_env, field_properties, flow_fields)
  end do

  if (is_root(par_env)) then
    print *, "Built ", size(flow_fields%fields), " dynamically-defined fields"
  end if

  call set_field_type(cell_centred_central, field_properties)
  call set_field_name("viscosity", field_properties)
  call create_field(par_env, field_properties, flow_fields) 
  call set_field_name("density", field_properties)
  call create_field(par_env, field_properties, flow_fields) 

  call set_vector_location(face, vec_properties)
  call set_size(par_env, mesh, vec_properties)
  call set_field_vector_properties(vec_properties, field_properties)
  call set_field_type(face_centred, field_properties)
  call set_field_name("mf", field_properties)
  call create_field(par_env, field_properties, flow_fields)

  call get_field(flow_fields, "u", u)
  call get_field(flow_fields, "v", v)
  call get_field(flow_fields, "mf", mf)
  call get_field(flow_fields, "viscosity", viscosity)
  call get_field(flow_fields, "density", density)
  call get_field(flow_fields, "scalar", scalar)

  ! Add fields to output list
  call add_field_to_outputlist(u)
  call add_field_to_outputlist(v)
  call add_field_to_outputlist(scalar)

  nullify(u)
  nullify(v)
  nullify(mf)
  nullify(viscosity)
  nullify(density)
  nullify(scalar)

  ! Initialise velocity field
  if (irank == par_env%root) print *, "Initialise velocity field"
  call initialise_flow(flow_fields)

  ! Initialise with default values
  if (irank == par_env%root) print *, "Initialise mat"
  call initialise(mat_properties)
  call initialise(scalar_equation_system)

  ! Create stiffness matrix
  call set_size(par_env, mesh, mat_properties)
  call set_nnz(5, mat_properties)
  call create_matrix(mat_properties, M)

  ! Create right-hand-side and solution vectors
  call initialise(vec_properties)
  call set_vector_location(cell, vec_properties)
  call set_size(par_env, mesh, vec_properties)
  call create_vector(vec_properties, source)

  call zero(source)
  call zero(M)

  call get_field(flow_fields, "scalar", scalar)
  call get_field(flow_fields, "mf", mf)
  call get_field(flow_fields, "density", density)
  call get_field(flow_fields, "viscosity", viscosity)
  call compute_fluxes(scalar, mf, viscosity, density, direction, M, source)
  nullify(scalar)
  nullify(mf)
  nullify(density)
  nullify(viscosity)

  call update(M) ! parallel assembly for M
  call update(source) ! parallel assembly for source

  call finalise(M)

  if(restart) then
    if (is_root(par_env)) then
      print*, "restart capability activated"
    end if
    call read_solution(par_env, case_path, mesh, flow_fields)
  end if 

  ! Create linear solver & set options
  if (irank == par_env%root) print *, "Solve"
  call get_field(flow_fields, "scalar", scalar)

  if(.not.unsteady) then
    num_steps = 1
    print*, "steady-state activated"
  else
    print*, "unsteady-state activated"
  end if

  do t = 1, num_steps
    call timer_register_start("Solver time inc I/O", timer_index_sol)
    call set_equation_system(par_env, source, scalar%values, M, scalar_equation_system)
    call create_solver(scalar_equation_system, scalar_solver)
    call solve(scalar_solver)
    nullify(scalar)

    if ((t == 1) .or. (t == num_steps) .or. (mod(t, write_frequency) == 0)) then
      if(.not. unsteady) then
        call write_solution(par_env, case_path, mesh, flow_fields)
      else
        call timer_start(timer_index_io_sol)
        !call write_mesh(par_env, case_path, mesh)
        call write_solution(par_env, case_path, mesh, flow_fields, t, num_steps, dt)
        call timer_stop(timer_index_io_sol)
      end if 
    end if 

    call timer_stop(timer_index_sol)

  end do

  ! Clean up
  deallocate(source)
  deallocate(M)
  deallocate(scalar_solver)

  call timer(end_time)

  if (par_env%proc_id == par_env%root) then
    print *, "Elapsed time = ", (end_time - start_time)
  end if

  call nullify_mesh_object()
  call cleanup_parallel_environment(par_env)
contains

  subroutine initialise_flow(flow_fields)
    use constants, only: insert_mode, ndim
    use types, only: vector_values, cell_locator, face_locator, neighbour_locator
    use meshing, only: create_cell_locator, get_global_index, count_neighbours, create_neighbour_locator, &
                       get_local_index, create_face_locator, get_local_index, get_face_normal, get_centre, &
                       get_local_num_cells
    use fv, only: calc_cell_coords
    use utils, only: clear_entries, set_mode, set_row, set_entry, set_values
    use vec, only: get_vector_data, restore_vector_data, create_vector_values

    ! Arguments
    type(fluid), intent(inout) :: flow_fields

    ! Local variables
    class(field), pointer :: u, v, mf, viscosity, density
    class(field), pointer :: scalar
    integer(ccs_int) :: n, count
    integer(ccs_int) :: n_local
    integer(ccs_int) :: index_p, global_index_p, index_f, index_nb
    real(ccs_real) :: u_val, v_val, scalar_val
    type(cell_locator) :: loc_p
    type(face_locator) :: loc_f
    type(neighbour_locator) :: loc_nb
    type(vector_values) :: u_vals, v_vals, scalar_vals
    real(ccs_real), dimension(:), pointer :: mf_data, viscosity_data, density_data

    real(ccs_real), dimension(ndim) :: x_p, x_f
    real(ccs_real), dimension(ndim) :: face_normal

    integer(ccs_int) :: nnb
    integer(ccs_int) :: j

    ! Get field pointers
    call get_field(flow_fields, "u", u)
    call get_field(flow_fields, "v", v)
    call get_field(flow_fields, "mf", mf)
    call get_field(flow_fields, "viscosity", viscosity)
    call get_field(flow_fields, "density", density)
    call get_field(flow_fields, "scalar", scalar)

    ! Set alias
    call get_local_num_cells(n_local)

    call create_vector_values(n_local, u_vals)
    call create_vector_values(n_local, v_vals)
    call create_vector_values(n_local, scalar_vals)
    call set_mode(insert_mode, u_vals)
    call set_mode(insert_mode, v_vals)
    call set_mode(insert_mode, scalar_vals)

    ! Set initial values for velocity fields
    do index_p = 1, n_local
      call create_cell_locator(index_p, loc_p)
      call get_global_index(loc_p, global_index_p)

      call get_centre(loc_p, x_p)

      u_val = 1.0_ccs_real !x_p(1) / domain_size
      v_val = 0.0_ccs_real !-x_p(2) / domain_size

      scalar_val = 0.0_ccs_real !x_p(1) 

      call set_row(global_index_p, u_vals)
      call set_entry(u_val, u_vals)
      call set_row(global_index_p, v_vals)
      call set_entry(v_val, v_vals)
      call set_row(global_index_p, scalar_vals)
      call set_entry(scalar_val, scalar_vals)
    end do

    call set_values(u_vals, u%values)
    call set_values(v_vals, v%values)
    call set_values(scalar_vals, scalar%values)

    deallocate (u_vals%global_indices)
    deallocate (v_vals%global_indices)
    deallocate (scalar_vals%global_indices)
    deallocate (u_vals%values)
    deallocate (v_vals%values)
    deallocate (scalar_vals%values)

    call get_vector_data(mf%values, mf_data)

    count = 0
    n = 0

    ! Loop over local cells and faces
    call get_local_num_cells(n_local)
    do index_p = 1, n_local

      call create_cell_locator(index_p, loc_p)
      call count_neighbours(loc_p, nnb)
      do j = 1, nnb

        call create_neighbour_locator(loc_p, j, loc_nb)
        call get_local_index(loc_nb, index_nb)

        ! if neighbour index is greater than previous face index
        if (index_nb > index_p) then ! XXX: abstract this test

          call create_face_locator(index_p, j, loc_f)
          call get_local_index(loc_f, index_f)
          call get_face_normal(loc_f, face_normal)
          call get_centre(loc_f, x_f)

          ! compute initial value based on current face coordinates
          ! mf_data(index_f) = 0.0_ccs_real * face_normal(1)
          mf_data(index_f) = 1.0_ccs_real * abs(face_normal(1))
        end if

      end do
    end do

    call restore_vector_data(mf%values, mf_data)

    call get_vector_data(viscosity%values, viscosity_data)
    viscosity_data(:) =  1.e-2_ccs_real
    call restore_vector_data(viscosity%values, viscosity_data)

    call get_vector_data(density%values, density_data)
    density_data(:) = 1.0_ccs_real
    call restore_vector_data(density%values, density_data)

    call update(u%values)
    call update(v%values)
    call update(scalar%values)
    call update(mf%values)
    call update(viscosity%values)
    call update(density%values)

    nullify(u)
    nullify(v)
    nullify(mf)
    nullify(viscosity)
    nullify(density)
    nullify(scalar)

  end subroutine initialise_flow

  ! Read YAML configuration file
  subroutine read_configuration(config_filename)

    use read_config, only: get_reference_number, get_value, &
                           get_relaxation_factors

    character(len=*), intent(in) :: config_filename

    class(*), pointer :: config_file  !< Pointer to CCS config file
    character(:), allocatable :: error

    config_file => parse(config_filename, error)
    if (allocated(error)) then
      call error_abort(trim(error))
    end if

    call get_variables(config_file, variable_names)
    if (size(variable_names) == 0) then
      call error_abort("No variables were specified.")
    end if
    print*,"no. of variables=",size(variable_names)
    call get_variable_types(config_file, variable_types)
    if (size(variable_types) /= size(variable_names)) then
       call error_abort("The number of variable types does not match the number of named variables")
    end if

    call get_value(config_file, 'restart', restart)

    call get_value(config_file, 'unsteady', unsteady)

    call get_value(config_file, 'iterations', num_iters)
    if (num_iters == huge(0)) then
      call error_abort("No value assigned to num_iters.")
    end if

    if(unsteady) then
      call get_value(config_file, 'steps', num_steps)
      if (num_steps == huge(0)) then
        call error_abort("No value assigned to num_steps.")
      end if

      call get_value(config_file, 'dt', dt)
      if (dt == huge(0.0)) then
        call error_abort("No value assigned to dt.")
      end if

      call get_value(config_file, 'write_frequency', write_frequency)
      if (write_frequency == huge(0.0)) then
        call error_abort("No value assigned to write_frequency.")
      end if
    end if 

    if (cps == huge(0)) then ! cps was not set on the command line
      call get_value(config_file, 'cps', cps)
      if (cps == huge(0)) then
        call error_abort("No value assigned to cps.")
      end if
    end if

    call get_value(config_file, 'L', domain_size)
    if (domain_size == huge(0.0)) then
      call error_abort("No value assigned to domain_size.")
    end if

    call get_value(config_file, 'target_residual', res_target)
    if (res_target == huge(0.0)) then
      call error_abort("No value assigned to target residual.")
    end if

    call get_relaxation_factors(config_file, u_relax=velocity_relax, p_relax=pressure_relax)
    if (velocity_relax == huge(0.0) .and. pressure_relax == huge(0.0)) then
      call error_abort("No values assigned to velocity and pressure underrelaxation.")
    end if

  end subroutine


end program scalar_advection
