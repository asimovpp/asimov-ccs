!> Program file for Sandia case
program sandia
#include "ccs_macros.inc"

  use petscvec
  use petscsys

  use core
  use ccs_base, only: mesh
  use case_config, only: num_steps, num_iters, dt, domain_size, write_frequency, &
                         velocity_relax, pressure_relax, res_target, case_name, &
                         write_gradients, velocity_solver_method_name, velocity_solver_precon_name, &
                         pressure_solver_method_name, pressure_solver_precon_name
  use constants, only: cell, face, ccsconfig, ccs_string_len, geoext, adiosconfig, ndim, &
                       cell_centred_central, cell_centred_upwind, face_centred, &
                       ccs_split_type_shared, ccs_split_type_low_high, ccs_split_undefined
  use meshing, only: set_mesh_object, nullify_mesh_object
  use kinds, only: ccs_real, ccs_int, ccs_long
  use types, only: field, field_spec, upwind_field, central_field, face_field, ccs_mesh, &
                   vector_spec, ccs_vector, io_environment, io_process, &
                   field_ptr, fluid, bc_profile
  use fields, only: create_field, set_field_config_file, set_field_n_boundaries, set_field_name, &
                    set_field_type, set_field_vector_properties, set_field_store_residuals, set_field_enable_cell_corrections
  use fortran_yaml_c_interface, only: parse
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments, sync, is_root
  use parallel_types, only: parallel_environment
  use vec, only: create_vector, set_vector_location
  use petsctypes, only: vector_petsc
  use scalars, only: update_scalars
  use utils, only: set_size, initialise, update, exit_print, &
                   add_field_to_outputlist, get_field, add_field, &
                   set_is_field_solved, &
                   allocate_fluid_fields
  use boundary_conditions, only: read_bc_config, allocate_bc_arrays, set_bc_profile
  use read_config, only: get_variables, get_case_name, &
                          get_variable_types
  use timestepping, only: set_timestep, activate_timestepping, initialise_old_values
  use mesh_utils, only: read_mesh, write_mesh
  use partitioning, only: compute_partitioner_input, &
                          partition_kway, compute_connectivity
  use fv, only: update_gradient
  use utils, only: str
  use timers, only: timer_init, timer_register_start, timer_register, timer_start, timer_stop, &
                    timer_print, timer_get_time, timer_print_all, timer_export_csv

  implicit none

  class(parallel_environment), allocatable:: par_env
  class(parallel_environment), allocatable:: shared_env
  character(len=:), allocatable:: input_path  ! Path to input directory
  character(len=:), allocatable:: case_path  ! Path to input directory with case name appended
  character(len = ccs_string_len), dimension(:), allocatable:: variable_names  ! variable names for BC reading
  integer(ccs_int), dimension(:), allocatable:: variable_types              ! cell centred upwind, central, etc.

  type(ccs_options) :: run_options

  integer(ccs_int):: it_start, it_end
  integer(ccs_int):: irank  ! MPI rank ID
  integer(ccs_int):: isize  ! Size of MPI world

  integer(ccs_int):: timer_index_total
  integer(ccs_int):: timer_index_init
  integer(ccs_int):: n_variables

  type(fluid):: flow_fields
  ! type(bc_profile), allocatable:: profile

  ! Launch MPI
  call initialise_parallel_environment(par_env)

  ! XXX: These values should be set by configuration (i.e. call after config)
  run_options%split_type = ccs_split_type_shared
  run_options%use_mpi_splitting = .true.
  call configure_parallelism(run_options, par_env, shared_env)
  irank = par_env%proc_id
  isize = par_env%num_procs

  call read_command_line_arguments(par_env, case_name = case_name, in_dir = input_path)

  ! XXX: set case name (should be absorbed into generic command line options reading later)
  run_options%init_mesh_type = build_mesh_3d
  run_options%case_name = case_name
  run_options%case_path = case_path
  run_options%mesh_path = case_name // "_mesh" // geoext
  run_options%cps = 10

  if (allocated(input_path)) then
    case_path = input_path // "/" // case_name
  else
    case_path = case_name
  end if

  run_options%config_file = case_path // ccsconfig

  call timer_register_start("Elapsed time", timer_index_total, is_total_time=.true.)

  call timer_register_start("Init time", timer_index_init)

  ! Read case name and runtime parameters from configuration file
  call read_configuration(run_options%config_file)

  run_options%variable_names = variable_names
  run_options%variable_types = variable_types

  if (is_root(par_env)) print *, "Starting ", case_name, " case!"

  ! set solver and preconditioner info
  velocity_solver_method_name = "gmres"
  velocity_solver_precon_name = "bjacobi"
  pressure_solver_method_name = "cg"
  pressure_solver_precon_name = "gamg"

  ! Set start and end iteration numbers (read from input file)
  it_start = 1
  it_end = num_iters

  ! Read mesh from .geo file
  call initialise_mesh(par_env, shared_env, run_options)

  ! Initialise fields
  if (irank == par_env%root) print *, "Initialise fields"

  ! Write gradients to solution file
  write_gradients = .false.

  ! Initialise the fields
  n_variables = size(run_options%variable_names)
  allocate(run_options%output(n_variables))
  allocate(run_options%solve(n_variables))
  run_options%output(1:4) = .true.
  run_options%output(5:7) = .false.
  run_options%output(8) = .true.
  run_options%solve(1:4) = .true.
  run_options%solve(5:) = .false.
  call initialise_fields(par_env, run_options, flow_fields)

  ! Initialise velocity field
  if (irank == par_env%root) print *, "Initialise velocity field"
  call initialise_flow(flow_fields, get_init_flow, get_init_mass_flux)

  ! Solve using SIMPLE algorithm
  if (irank == par_env%root) print *, "Start SIMPLE"

  ! Print the run configuration
  if (irank == par_env%root) then
    call print_configuration()
  end if

  call activate_timestepping()
  call set_timestep(dt)

  call timer_stop(timer_index_init)

  ! XXX: These values should be set during configuration
  run_options%case_path = case_path
  run_options%num_steps = num_steps
  run_options%dt = dt
  run_options%write_frequency = write_frequency
  run_options%it_start = it_start
  run_options%it_end = it_end
  run_options%res_target = res_target
  call run_solver(par_env, run_options, flow_fields)
  
  ! Clean-up

  call timer_stop(timer_index_total)

  call timer_print_all(par_env)
  call timer_export_csv(par_env)

  ! Finalise MPI
  call cleanup_parallel_environment(par_env)

contains

  ! Read YAML configuration file
  subroutine read_configuration(config_filename)

    use read_config, only: get_reference_number, get_value, &
                           get_relaxation_factors

    character(len=*), intent(in):: config_filename

    class(*), pointer:: config_file  !< Pointer to CCS config file
    character(:), allocatable:: error

    config_file => parse(config_filename, error)
    if (allocated(error)) then
      call error_abort(trim(error))
    end if

    call get_variables(config_file, variable_names)
    if (size(variable_names) == 0) then
      call error_abort("No variables were specified.")
    end if
    call get_variable_types(config_file, variable_types)
    if (size(variable_types) /= size(variable_names)) then
       call error_abort("The number of variable types does not match the number of named variables")
    end if


    call get_value(config_file, 'steps', num_steps)
    if (num_steps == huge(0)) then
      call error_abort("No value assigned to num_steps.")
    end if

    call get_value(config_file, 'iterations', num_iters)
    if (num_iters == huge(0)) then
      call error_abort("No value assigned to num_iters.")
    end if

    call get_value(config_file, 'dt', dt)
    if (dt == huge(0.0)) then
      call error_abort("No value assigned to dt.")
    end if

    call get_value(config_file, 'write_frequency', write_frequency)
    if (write_frequency == huge(0.0)) then
      call error_abort("No value assigned to write_frequency.")
    end if

    call get_value(config_file, 'L', domain_size)
    if (domain_size == huge(0.0)) then
      call error_abort("No value assigned to domain_size.")
    end if

    call get_value(config_file, 'target_residual', res_target)
    if (res_target == huge(0.0)) then
      call error_abort("No value assigned to target residual.")
    end if

    call get_relaxation_factors(config_file, u_relax = velocity_relax, p_relax = pressure_relax)
    if (velocity_relax == huge(0.0) .and. pressure_relax == huge(0.0)) then
      call error_abort("No values assigned to velocity and pressure underrelaxation.")
    end if

  end subroutine

  ! Print test case configuration
  subroutine print_configuration()

    ! XXX: this should eventually be replaced by something nicely formatted that uses "write"
    print *, " "
    print *, "******************************************************************************"
    print *, "* Solving the ", case_name, " case"
    print *, "******************************************************************************"
    print *, " "
    print *, "******************************************************************************"
    print *, "* SIMULATION LENGTH"
    print *, "* Running for ", num_steps, "timesteps and ", num_iters, "iterations"
    write (*, '(1x, a, e10.3)') "* Time step size: ", dt
    print *, "******************************************************************************"
    print *, "* MESH SIZE"
    print *, "* Global number of cells is ", mesh%topo%global_num_cells
    print *, "******************************************************************************"
    print *, "* RELAXATION FACTORS"
    write (*, '(1x, a, e10.3)') "* velocity: ", velocity_relax
    write (*, '(1x, a, e10.3)') "* pressure: ", pressure_relax
    print *, "******************************************************************************"

  end subroutine


  pure subroutine get_init_flow(loc_p, field_name, init_val)
    use types, only: cell_locator
    use meshing, only: get_centre
    type(cell_locator), intent(in) :: loc_p
    character(len=*), intent(in) :: field_name
    real(ccs_real), intent(inout) :: init_val
    real(ccs_real), dimension(ndim) :: x_p

    if (field_name == "scalar") then
      call get_centre(loc_p, x_p)
      if (x_p(1) < -0.08) then
        init_val = 1.0_ccs_real 
      else
        init_val = 0.0_ccs_real 
      end if
    else ! anything but scalar field
      init_val = 0.0_ccs_real
    end if

  end subroutine
  
  pure subroutine get_init_mass_flux(loc_f, init_val)
    use types, only: face_locator
    type(face_locator), intent(in) :: loc_f
    real(ccs_real), intent(inout) :: init_val

    associate (foo => loc_f, bar => init_val)

    end associate

  end subroutine

end program sandia
