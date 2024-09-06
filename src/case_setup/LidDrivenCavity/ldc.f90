!v Program file for LidDrivenCavity case
!
!  @build mpi+petsc

program ldc
#include "ccs_macros.inc"

  use petscvec
  use petscsys

  use ccs_base, only: mesh
  use case_config, only: num_steps, num_iters, dt, cps, domain_size, case_name, &
                         velocity_relax, pressure_relax, res_target, &
                         write_gradients, velocity_solver_method_name, velocity_solver_precon_name, &
                         pressure_solver_method_name, pressure_solver_precon_name, restart, unsteady, write_frequency
  use constants, only: cell, face, ccsconfig, ccs_string_len, &
                       cell_centred_central, cell_centred_upwind, face_centred, &
                       ccs_split_type_shared, ccs_split_type_low_high
  use kinds, only: ccs_real, ccs_int
  use types, only: field, field_spec, upwind_field, central_field, gamma_field, face_field, ccs_mesh, &
                   vector_spec, ccs_vector, field_ptr, fluid
  use fields, only: create_field, set_field_config_file, set_field_n_boundaries, set_field_name, &
                    set_field_type, set_field_vector_properties, set_field_store_residuals
  use fortran_yaml_c_interface, only: parse
  use parallel, only: initialise_parallel_environment, &
                      create_new_par_env, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments, sync, is_root
  use meshing, only: set_mesh_object, nullify_mesh_object
  use parallel_types, only: parallel_environment
  use mesh_utils, only: build_mesh, write_mesh, build_square_mesh
  use meshing, only: get_global_num_cells, get_local_num_cells
  use vec, only: create_vector, set_vector_location, get_vector_data, restore_vector_data
  use petsctypes, only: vector_petsc
  use pv_coupling, only: solve_nonlinear
  use utils, only: set_size, initialise, update, exit_print, add_field_to_outputlist, &
                   get_field, set_is_field_solved, &
                   allocate_fluid_fields, dealloc_fluid_fields
  use boundary_conditions, only: read_bc_config, allocate_bc_arrays
  use read_config, only: get_variables, get_boundary_count, get_boundary_names, get_store_residuals, get_variable_types
  use io_visualisation, only: write_solution, read_solution
  use timers, only: timer_init, timer_register_start, timer_register, timer_start, timer_stop, timer_print, timer_print_all

  implicit none

  class(parallel_environment), allocatable :: par_env
  class(parallel_environment), allocatable :: shared_env
  character(len=:), allocatable :: input_path  ! Path to input directory
  character(len=:), allocatable :: case_path  ! Path to input directory with case name appended
  character(len=:), allocatable :: ccs_config_file ! Config file for CCS
  character(len=ccs_string_len), dimension(:), allocatable :: variable_names  ! variable names for BC reading
  integer(ccs_int), dimension(:), allocatable :: variable_types              ! cell centred upwind, central, etc.

  type(vector_spec) :: vec_properties

  type(field_spec) :: field_properties
  class(field), pointer :: u, v, w, p, mf, viscosity, density

  integer(ccs_int) :: n_boundaries

  integer(ccs_int) :: it_start, it_end
  integer(ccs_int) :: irank ! MPI rank ID
  integer(ccs_int) :: isize ! Size of MPI world

  integer(ccs_int) :: timer_index_init, timer_index_total, timer_index_sol
  integer(ccs_int) :: i

  logical :: u_sol = .true.  ! Default equations to solve for LDC case
  logical :: v_sol = .true.
  logical :: w_sol = .true.
  logical :: p_sol = .true.

  logical :: store_residuals

  type(fluid) :: flow_fields

  logical :: use_mpi_splitting

  character(len=128), dimension(:), allocatable :: bnd_names

  integer(ccs_int):: t          ! Timestep counter
  integer(ccs_int):: timer_index_io_sol
  
  ! Launch MPI
  call initialise_parallel_environment(par_env)
  call timer_init()

  irank = par_env%proc_id
  isize = par_env%num_procs

  ! Create shared memory communicator for each node
  use_mpi_splitting = .true.
  call create_new_par_env(par_env, ccs_split_type_shared, use_mpi_splitting, shared_env)

  call read_command_line_arguments(par_env, cps, case_name=case_name, in_dir=input_path)

  if (allocated(input_path)) then
    case_path = input_path // "/" // case_name
  else
    case_path = case_name
  end if

  ccs_config_file = case_path // ccsconfig

  call timer_register_start("Elapsed time", timer_index_total, is_total_time=.true.)
  call timer_register_start("Init time", timer_index_init)

  ! Read case name from configuration file
  call read_configuration(ccs_config_file)

  if (irank == par_env%root) print *, "Starting ", case_name, " case!"

  ! set solver and preconditioner info
  velocity_solver_method_name = "gmres"
  velocity_solver_precon_name = "bjacobi"
  pressure_solver_method_name = "cg"
  pressure_solver_precon_name = "gamg"

  ! Set start and end iteration numbers (read from input file)
  it_start = 1
  it_end = num_iters
  
  ! Create a mesh
  if (irank == par_env%root) print *, "Building mesh"
  call get_boundary_names(ccs_config_file, bnd_names)
  !mesh = build_mesh(par_env, cps, cps, cps, 1.0_ccs_real, bnd_names)   ! 3-D mesh
  mesh = build_square_mesh(par_env, shared_env, cps, 1.0_ccs_real, bnd_names)      ! 2-D mesh
  call set_mesh_object(mesh)

  ! Initialise fields
  if (irank == par_env%root) print *, "Initialise fields"

  ! Write gradients to solution file
  write_gradients = .false.

  ! Create and initialise field vectors
  call initialise(vec_properties)
  call get_boundary_count(ccs_config_file, n_boundaries)
  call get_store_residuals(ccs_config_file, store_residuals)

  call set_vector_location(cell, vec_properties)
  call set_size(par_env, mesh, vec_properties)
  call set_field_config_file(ccs_config_file, field_properties)
  call set_field_n_boundaries(n_boundaries, field_properties)
  call set_field_store_residuals(store_residuals, field_properties)

  call set_field_vector_properties(vec_properties, field_properties)

  if (is_root(par_env)) then
    print *, "Build field list"
  end if

  ! We expect to find u,v,w,p,p_prime
  do i = 1, size(variable_names)
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

  ! Get field pointers
  call get_field(flow_fields, "u", u)
  call get_field(flow_fields, "v", v)
  call get_field(flow_fields, "w", w)
  call get_field(flow_fields, "p", p)
  call get_field(flow_fields, "mf", mf)
  call get_field(flow_fields, "viscosity", viscosity)
  call get_field(flow_fields, "density", density)
  
  ! Add fields to output list
  call add_field_to_outputlist(u)
  call add_field_to_outputlist(v)
  call add_field_to_outputlist(w)
  call add_field_to_outputlist(p)

  ! Initialise velocity field
  if (irank == par_env%root) print *, "Initialise velocity field"
  call initialise_velocity(flow_fields) 

  ! XXX: This should get incorporated as part of create_field subroutines
  call set_is_field_solved(u_sol, u)
  call set_is_field_solved(v_sol, v)
  call set_is_field_solved(w_sol, w)
  call set_is_field_solved(p_sol, p)

  ! Nullify pointers for safety
  nullify(u)
  nullify(v)
  nullify(w)
  nullify(p)
  nullify(mf)
  nullify(viscosity)
  nullify(density)

  if(restart) then
    if (is_root(par_env)) then
      print*, "restart capability activated"
    end if
    call read_solution(par_env, case_path, mesh, flow_fields)
  end if 

  if (irank == par_env%root) then
    call print_configuration()
  end if

  call timer_stop(timer_index_init)
  call timer_register("I/O time for solution", timer_index_io_sol)

  if(.not.unsteady) then
    num_steps = 1
    print*, "steady-state activated"
  else
    print*, "unsteady-state activated"
  end if

  do t = 1, num_steps
    call timer_register_start("Solver time inc I/O", timer_index_sol)
    ! Solve using SIMPLE algorithm
    if (irank == par_env%root) print *, "Start SIMPLE"
    call solve_nonlinear(par_env, mesh, eval_sources, it_start, it_end, res_target, &
                     flow_fields)

    ! Write out mesh and solution
    call write_mesh(par_env, case_path, mesh)
    if (par_env%proc_id == par_env%root) then
      print *, "TIME = ", t
    end if

    if ((t == 1) .or. (t == num_steps) .or. (mod(t, write_frequency) == 0)) then
      if(.not. unsteady) then
        call write_solution(par_env, case_path, mesh, flow_fields)
      else 
        call timer_start(timer_index_io_sol)
        call write_solution(par_env, case_path, mesh, flow_fields, t, num_steps, dt)
        call timer_stop(timer_index_io_sol)
      end if 
    end if

    call timer_stop(timer_index_sol)
  end do

  ! Clean-up
  call dealloc_fluid_fields(flow_fields)

  call timer_stop(timer_index_total)

  call timer_print_all(par_env)

  call nullify_mesh_object()

  ! Finalise MPI
  call cleanup_parallel_environment(par_env)

contains

  ! Read YAML configuration file
  subroutine read_configuration(config_filename)

    use read_config, only: get_value, get_relaxation_factors

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

    call get_value(config_file, 'iterations', num_iters) ! steady-state
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

  ! Print test case configuration
  subroutine print_configuration()

    integer(ccs_int) :: global_num_cells

    call get_global_num_cells(global_num_cells)

    ! XXX: this should eventually be replaced by something nicely formatted that uses "write"
    print *, " "
    print *, "******************************************************************************"
    print *, "* Solving the ", case_name, " case"
    print *, "******************************************************************************"
    print *, " "
    print *, "******************************************************************************"
    print *, "* SIMULATION LENGTH"
    if (unsteady) then
      print *, "* Running for ", num_steps, "timesteps and ", num_iters, "iterations"
      write (*, '(1x, a, e10.3)') "* Time step size: ", dt
    else
      print *, "* Running for ", num_iters, "iterations"
    end if 

    print *, "******************************************************************************"
    print *, "* MESH SIZE"
    print *, "* Cells per side: ", cps
    write (*, '(1x,a,e10.3)') "* Domain size: ", domain_size
    print *, "* Global number of cells is ", global_num_cells
    print *, "******************************************************************************"
    print *, "* RELAXATION FACTORS"
    write (*, '(1x,a,e10.3)') "* velocity: ", velocity_relax
    write (*, '(1x,a,e10.3)') "* pressure: ", pressure_relax
    print *, "******************************************************************************"

  end subroutine

  subroutine initialise_velocity(flow_fields)

    use constants, only: add_mode
    use types, only: vector_values, cell_locator
    use meshing, only: create_cell_locator, get_global_index, get_local_num_cells
    use fv, only: calc_cell_coords
    use utils, only: clear_entries, set_mode, set_row, set_entry, set_values
    use vec, only: get_vector_data, restore_vector_data, create_vector_values

    ! Arguments
    type(fluid), intent(inout) :: flow_fields

    ! Local variables
    class(field), pointer :: u, v, w, mf, mu, rho
    integer(ccs_int) :: row, col
    integer(ccs_int) :: index_p, global_index_p, n_local
    real(ccs_real) :: u_val, v_val, w_val
    type(cell_locator) :: loc_p
    type(vector_values) :: u_vals, v_vals, w_vals
    real(ccs_real), dimension(:), pointer :: mf_data, viscosity_data, density_data
    
    ! Set alias
    call get_local_num_cells(n_local)
    call create_vector_values(n_local, u_vals)
    call create_vector_values(n_local, v_vals)
    call create_vector_values(n_local, w_vals)
    call set_mode(add_mode, u_vals)
    call set_mode(add_mode, v_vals)
    call set_mode(add_mode, w_vals)

    ! Set initial values for velocity fields
    do index_p = 1, n_local
      call create_cell_locator(index_p, loc_p)
      call get_global_index(loc_p, global_index_p)
      call calc_cell_coords(global_index_p, cps, row, col)

      u_val = 0.0_ccs_real
      v_val = 0.0_ccs_real
      w_val = 0.0_ccs_real

      call set_row(global_index_p, u_vals)
      call set_entry(u_val, u_vals)
      call set_row(global_index_p, v_vals)
      call set_entry(v_val, v_vals)
      call set_row(global_index_p, w_vals)
      call set_entry(w_val, w_vals)
    end do

    call get_field(flow_fields, "u", u)
    call get_field(flow_fields, "v", v)
    call get_field(flow_fields, "w", w)
    
    call set_values(u_vals, u%values)
    call set_values(v_vals, v%values)
    call set_values(w_vals, w%values)

    call update(u%values)
    call update(v%values)
    call update(w%values)
    nullify(u)
    nullify(v)
    nullify(w)

    deallocate (u_vals%global_indices)
    deallocate (v_vals%global_indices)
    deallocate (w_vals%global_indices)
    deallocate (u_vals%values)
    deallocate (v_vals%values)
    deallocate (w_vals%values)

    call get_field(flow_fields, "mf", mf)
    call get_field(flow_fields, "viscosity", mu)
    call get_field(flow_fields, "density", rho)

    call get_vector_data(mf%values, mf_data)
    mf_data(:) = 0.0_ccs_real
    call restore_vector_data(mf%values, mf_data)

    call get_vector_data(mu%values, viscosity_data)
    viscosity_data(:) =  1.e-2_ccs_real
    call restore_vector_data(mu%values, viscosity_data)

    call get_vector_data(rho%values, density_data)
    density_data(:) = 1.0_ccs_real
    call restore_vector_data(rho%values, density_data)

    call update(mf%values)
    call update(mu%values)
    call update(rho%values)
    nullify(mf)
    nullify(mu)
    nullify(rho)

  end subroutine initialise_velocity

  !> Case-specific source terms
  subroutine eval_sources(flow, phi, R, S)
    use types, only: fluid, field, ccs_vector
    use fv, only: zero_sources

    type(fluid), intent(in) :: flow !< Provides access to full flow field
    class(field), intent(in) :: phi !< Field being transported
    class(ccs_vector), intent(inout) :: R !< Work vector (for evaluating linear/implicit sources)
    class(ccs_vector), intent(inout) :: S !< Work vector (for evaluating fixed/explicit sources)
    
    ! Dummy implementation - just zeros the sources, see sero_sources for example implementation
    call zero_sources(flow, phi, R, S)
    
  end subroutine eval_sources

end program ldc
