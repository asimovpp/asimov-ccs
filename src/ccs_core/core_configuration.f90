submodule(core) core_configuration
#include "ccs_macros.inc"

  use read_config, only: get_variable_types, &
                         get_value, &
                         get_variables, get_relaxation_factors
  use utils, only: exit_print

  implicit none

contains

  !v Subroutine to get the runtime configuration.
  module subroutine get_config(par_env, run_options)

    use constants, only: ccs_split_type_shared
    use parallel, only: read_command_line_arguments

    class(parallel_environment), intent(in) :: par_env
    type(ccs_options), intent(out) :: run_options

    character(len=:), allocatable :: case_name
    character(len=:), allocatable :: input_path
    
    call read_command_line_arguments(par_env, case_name = case_name, in_dir = input_path)
    call set_case_paths(input_path, case_name, run_options)   

    ! Read case name and runtime parameters from configuration file
    call read_configuration(run_options)

    ! XXX: These values should be set by configuration (i.e. call after config)
    run_options%split_type = ccs_split_type_shared
    run_options%use_mpi_splitting = .true.

    ! XXX: These values should be set during configuration

    ! Print the run configuration
    call print_configuration(par_env, run_options)

  end subroutine get_config

  ! Set the path names used by the CCS case
  subroutine set_case_paths(input_path, case_name, run_options)

    use constants, only: ccsconfig
    
    character(len=:), allocatable, intent(in) :: input_path
    character(len=*), intent(in) :: case_name
    type(ccs_options), intent(inout) :: run_options

    character(len=:), allocatable :: case_path
    character(len=*), parameter :: geoext = ".geo"
    
    run_options%paths%case_name = case_name
    run_options%paths%mesh_path = case_name // "_mesh" // geoext
    if (allocated(input_path)) then
      case_path = input_path // "/" // case_name
    else
      case_path = case_name
    end if
    run_options%paths%case_path = case_path

    run_options%paths%ccs_config_file = case_path // ccsconfig
    
  end subroutine set_case_paths
  
  ! Read YAML configuration file
  subroutine read_configuration(run_options)

    use fortran_yaml_c_interface, only: parse

    type(ccs_options), intent(inout) :: run_options

    class(*), pointer:: config_file  !< Pointer to CCS config file
    character(:), allocatable:: error

    config_file => parse(run_options%paths%ccs_config_file, error)
    if (allocated(error)) then
      call error_abort(trim(error))
    end if

    call get_variable_definitions(config_file, run_options)
    call get_solver_options(config_file, run_options)
    call get_io_options(config_file, run_options)
    call get_mesh_options(config_file, run_options)

  end subroutine

  subroutine get_mesh_options(config_file, run_options)
    
    class(*), pointer, intent(in) :: config_file
    type(ccs_options), intent(inout) :: run_options

    real(ccs_real) :: domain_size

    run_options%mesh%init_mesh_type = build_mesh_3d

    if (run_options%mesh%init_mesh_type == read_input_mesh) then
      call error_abort("Reading the mesh hasn't been implemented yet!")
    else
      run_options%mesh%cps = 10

      domain_size = huge(0.0_ccs_real)
      call get_value(config_file, 'L', domain_size)
      if (domain_size == huge(0.0_ccs_real)) then
        call error_abort("No value assigned to domain_size.")
      end if
      run_options%mesh%domain_size = domain_size
    end if
    
  end subroutine get_mesh_options

  subroutine get_io_options(config_file, run_options)
    
    class(*), pointer, intent(in) :: config_file
    type(ccs_options), intent(inout) :: run_options

    integer(ccs_int) :: write_frequency

    write_frequency = huge(0)
    call get_value(config_file, 'write_frequency', write_frequency)
    if (write_frequency == huge(0.0)) then
      call error_abort("No value assigned to write_frequency.")
    end if
    run_options%io%write_frequency = write_frequency

  end subroutine get_io_options

  subroutine get_variable_definitions(config_file, run_options)

    use constants, only: ccs_string_len
    
    class(*), pointer, intent(in) :: config_file
    type(ccs_options), intent(inout) :: run_options

    character(len=ccs_string_len), dimension(:), allocatable :: variable_names
    integer(ccs_int), dimension(:), allocatable :: variable_types

    allocate(variable_names(0)) ! Default size to 0 for error checking
    call get_variables(config_file, variable_names)
    if (size(variable_names) == 0) then
      call error_abort("No variables were specified.")
    end if
    call get_variable_types(config_file, variable_types)
    if (size(variable_types) /= size(variable_names)) then
       call error_abort("The number of variable types does not match the number of named variables")
    end if

    run_options%variables%variable_names = variable_names
    run_options%variables%variable_types = variable_types
    
  end subroutine get_variable_definitions

  subroutine get_solver_options(config_file, run_options)

    class(*), pointer, intent(in) :: config_file
    type(ccs_options), intent(inout) :: run_options

    integer(ccs_int) :: num_steps
    integer(ccs_int) :: num_iters
    real(ccs_real) :: dt
    real(ccs_real) :: res_target
    real(ccs_real) :: velocity_relax, pressure_relax

    ! Default values for error checking
    num_steps = huge(0)
    num_iters = huge(0)
    dt = huge(0.0_ccs_real)
    res_target = huge(0.0_ccs_real)
    velocity_relax = huge(0.0_ccs_real)
    pressure_relax = huge(0.0_ccs_real)
    
    call get_value(config_file, 'steps', num_steps)
    if (num_steps == huge(0)) then
      call error_abort("No value assigned to num_steps.")
    end if
    run_options%solve%num_steps = num_steps

    call get_value(config_file, 'iterations', num_iters)
    if (num_iters == huge(0)) then
      call error_abort("No value assigned to num_iters.")
    end if
    run_options%solve%num_iters = num_iters

    call get_value(config_file, 'dt', dt)
    if (dt == huge(0.0)) then
      call error_abort("No value assigned to dt.")
    end if
    run_options%solve%dt = dt

    call get_value(config_file, 'target_residual', res_target)
    if (res_target == huge(0.0)) then
      call error_abort("No value assigned to target residual.")
    end if
    run_options%solve%res_target = res_target

    call get_relaxation_factors(config_file, u_relax = velocity_relax, p_relax = pressure_relax)
    if (velocity_relax == huge(0.0) .and. pressure_relax == huge(0.0)) then
      call error_abort("No values assigned to velocity and pressure underrelaxation.")
    end if
    run_options%solve%velocity_relax = velocity_relax
    run_options%solve%pressure_relax = pressure_relax

    ! XXX: Are these still relevant?
    run_options%solve%it_start = 1
    run_options%solve%it_end = num_iters
    
  end subroutine get_solver_options

  ! Print test case configuration
  subroutine print_configuration(par_env, run_options)

    use ccs_base, only: mesh
    use parallel, only: is_root

    class(parallel_environment), intent(in) :: par_env
    type(ccs_options), intent(in) :: run_options

    if (is_root(par_env)) then
      associate(case_name => run_options%paths%case_name, &
           num_steps => run_options%solve%num_steps, &
           num_iters => run_options%solve%num_iters, &
           dt => run_options%solve%dt, &
           velocity_relax => run_options%solve%velocity_relax, &
           pressure_relax => run_options%solve%pressure_relax)
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
      end associate
    end if
    
  end subroutine
end submodule core_configuration
