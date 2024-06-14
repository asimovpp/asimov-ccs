submodule(core) core_configuration
#include "ccs_macros.inc"

  use read_config, only: get_variable_types, &
                         get_value, &
                         get_variables, get_relaxation_factors, &
                         get_output_type, get_solve, &
                         get_reference_number
  use utils, only: exit_print

  implicit none

  integer(ccs_int), save :: cps_cmdline = huge(0) ! Cells-per side, potentially read from commandline
contains

  !v Subroutine to get the runtime configuration.
  module subroutine get_config(par_env, run_options)

    use constants, only: ccs_split_type_shared
    use parallel, only: read_command_line_arguments

    class(parallel_environment), intent(in) :: par_env
    type(ccs_options), intent(out) :: run_options

    character(len=:), allocatable :: case_name
    character(len=:), allocatable :: input_path
    
    print *, "READ CMDLINE"
    call read_command_line_arguments(par_env, cps = cps_cmdline, case_name = case_name, &
                                     in_dir = input_path)
    print *, "SET PATHS"
    call set_case_paths(input_path, case_name, run_options%paths)

    ! XXX: These values should be set by configuration (i.e. call after config)
    run_options%parallel%split_type = ccs_split_type_shared
    run_options%parallel%use_mpi_splitting = .true.

    ! Read case name and runtime parameters from configuration file
    print *, "BUILD CONFIG"
    call read_configuration(run_options)

    ! Print the run configuration
    call print_configuration(par_env, run_options)

  end subroutine get_config

  ! Set the path names used by the CCS case
  subroutine set_case_paths(input_path, case_name, paths)

    use constants, only: ccsconfig
    
    character(len=:), allocatable, intent(in) :: input_path
    character(len=*), intent(in) :: case_name
    type(ccs_paths), intent(inout) :: paths

    character(len=:), allocatable :: case_path
    character(len=*), parameter :: geoext = ".geo"
    
    paths%case_name = case_name
    paths%mesh_path = case_name // "_mesh" // geoext
    if (allocated(input_path)) then
      case_path = input_path // "/" // case_name
    else
      case_path = case_name
    end if
    paths%case_path = case_path

    paths%ccs_config_file = case_path // ccsconfig
    
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

    print *, "+ VARS"
    call get_variable_definitions(config_file, run_options%variables)
    print *, "+ REFS"
    call get_reference_values(config_file, run_options%reference_values)
    print *, "+ SOL"
    call get_solver_options(config_file, run_options%solve)
    print *, "+ IO"
    call get_io_options(config_file, run_options%io)
    print *, "+ MSH"
    call get_mesh_options(config_file, run_options%mesh)
    print *, "+ DONE"

  end subroutine

  subroutine get_mesh_options(config_file, mesh_opt)
    
    class(*), pointer, intent(in) :: config_file
    type(mesh_options), intent(out) :: mesh_opt

    real(ccs_real) :: domain_size

    mesh_opt%init_mesh_type = build_mesh_2d
    if (mesh_opt%init_mesh_type == read_input_mesh) then
      call error_abort("Reading the mesh hasn't been implemented yet!")
    else
      if (cps_cmdline == huge(0)) then  ! cps was not set on the command line
        call get_value(config_file, 'cps', mesh_opt%cps)
        if (mesh_opt%cps == huge(0)) then
          call error_abort("No value assigned to cps.")
        end if
      else
        mesh_opt%cps = cps_cmdline
      end if

      domain_size = huge(0.0_ccs_real)
      call get_value(config_file, 'L', domain_size)
      if (domain_size == huge(0.0_ccs_real)) then
        call error_abort("No value assigned to domain_size.")
      end if
      mesh_opt%domain_size = domain_size
    end if

    call get_value(config_file, 'compute_bwidth', mesh_opt%compute_bwidth)
    call get_value(config_file, 'compute_partqual', mesh_opt%compute_partqual)
    
  end subroutine get_mesh_options

  subroutine get_io_options(config_file, io)
    
    class(*), pointer, intent(in) :: config_file
    type(io_options), intent(out) :: io

    integer(ccs_int) :: write_frequency

    write_frequency = huge(0)
    call get_value(config_file, 'write_frequency', write_frequency)
    if (write_frequency == huge(0.0)) then
      call error_abort("No value assigned to write_frequency.")
    end if
    io%write_frequency = write_frequency

  end subroutine get_io_options

  subroutine get_variable_definitions(config_file, variables)

    use constants, only: ccs_string_len
    
    class(*), pointer, intent(in) :: config_file
    type(variable_options), intent(out) :: variables

    character(len=ccs_string_len), dimension(:), allocatable :: variable_names
    integer(ccs_int), dimension(:), allocatable :: variable_types

    character(len=:), allocatable :: post_type ! Where is I/O (cell/vertex)? Currently unused

    ! allocate(variable_names(0)) ! Default size to 0 for error checking
    call get_variables(config_file, variable_names)
    if (size(variable_names) == 0) then
      call error_abort("No variables were specified.")
    end if
    call get_variable_types(config_file, variable_types)
    if (size(variable_types) /= size(variable_names)) then
       call error_abort("The number of variable types does not match the number of named variables")
    end if

    variables%variable_names = variable_names
    variables%variable_types = variable_types
    
    call get_output_type(config_file, post_type, variables%output_variables)
    call get_solve(config_file, variables%solved_variables)

  end subroutine get_variable_definitions

  subroutine get_reference_values(config_file, reference_values)

    class(*), pointer, intent(in) :: config_file
    type(ref_vals), intent(out) :: reference_values

    ! Set defaults
    reference_values%p_ref = 0.0
    reference_values%p_total = 1.01325e5_ccs_real ! 1 Atmosphere
    reference_values%temp_ref = 293.15_ccs_real          ! 20C (in Kelvin)
    reference_values%dens_ref = 1.19_ccs_real             ! Air @ STP
    reference_values%visc_ref = 1.0e-5_ccs_real         ! Air @ STP
    reference_values%velo_ref = 1.0_ccs_real
    reference_values%len_ref = 1.0_ccs_real
    reference_values%pref_at_cell = -1 ! To be ignored
    
    ! Read
    call get_reference_number(config_file, &
         reference_values%p_ref, reference_values%p_total, reference_values%temp_ref, &
         reference_values%dens_ref, reference_values%visc_ref, &
         reference_values%velo_ref, reference_values%len_ref, &
         reference_values%pref_at_cell)
    
  end subroutine get_reference_values
  
  subroutine get_solver_options(config_file, solve)

    class(*), pointer, intent(in) :: config_file
    type(solver_options), intent(out) :: solve

    integer(ccs_int) :: num_steps
    integer(ccs_int) :: num_iters
    real(ccs_real) :: dt
    real(ccs_real) :: res_target
    real(ccs_real) :: velocity_relax, pressure_relax
    logical :: present, required

    ! Default values for error checking
    num_steps = huge(0)
    num_iters = huge(0)
    dt = huge(0.0_ccs_real)
    res_target = huge(0.0_ccs_real)
    velocity_relax = huge(0.0_ccs_real)
    pressure_relax = huge(0.0_ccs_real)

    call get_value(config_file, 'iterations', num_iters)
    if (num_iters == huge(0)) then
      call error_abort("No value assigned to num_iters.")
    end if
    solve%num_iters = num_iters

    required = .false.
    call get_value(config_file, 'steps', num_steps, present, required)
    call get_value(config_file, 'dt', dt, present, required)
    if ((num_steps == huge(0)) .and. (dt == huge(0.0))) then
      print *, "Steady-state solver"
    else if ((num_steps == huge(0)) .or. (dt == huge(0.0))) then
      call error_abort("No value assigned to either num_steps or dt.")
    end if
    solve%num_steps = num_steps
    solve%dt = dt

    call get_value(config_file, 'target_residual', res_target)
    if (res_target == huge(0.0)) then
      call error_abort("No value assigned to target residual.")
    end if
    solve%res_target = res_target

    call get_relaxation_factors(config_file, u_relax = velocity_relax, p_relax = pressure_relax)
    if (velocity_relax == huge(0.0) .and. pressure_relax == huge(0.0)) then
      call error_abort("No values assigned to velocity and pressure underrelaxation.")
    end if
    solve%velocity_relax = velocity_relax
    solve%pressure_relax = pressure_relax

    ! XXX: Are these still relevant?
    solve%it_start = 1
    solve%it_end = num_iters

  end subroutine get_solver_options

  ! Print test case configuration
  subroutine print_configuration(par_env, run_options)

    use ccs_base, only: mesh
    use parallel, only: is_root

    class(parallel_environment), intent(in) :: par_env
    type(ccs_options), intent(in) :: run_options

    integer :: i

    if (is_root(par_env)) then
      associate(case_name => run_options%paths%case_name, &
           num_steps => run_options%solve%num_steps, &
           num_iters => run_options%solve%num_iters, &
           dt => run_options%solve%dt, &
           cps => run_options%mesh%cps, &
           domain_size => run_options%mesh%domain_size, &
           velocity_relax => run_options%solve%velocity_relax, &
           pressure_relax => run_options%solve%pressure_relax)
        ! XXX: this should eventually be replaced by something nicely formatted that uses "write"
        print *, " "
        print *, "******************************************************************************"
        print *, "* Solving the ", case_name, " case"
        print *, "******************************************************************************"
        print *, "Solved variables: "
        do i = 1, size(run_options%variables%solved_variables)
          print *, "- ", run_options%variables%solved_variables(i)
        end do
        print *, "******************************************************************************"
        print *, "* SIMULATION LENGTH"
        if (dt /= huge(dt)) then
          print *, "* Running for ", num_steps, "timesteps and ", num_iters, "iterations"
          write (*, '(1x, a, e10.3)') "* Time step size: ", dt
        else
          print *, "* Running for ", num_iters, "iterations"
        end if
        print *, "******************************************************************************"
        print *, "* MESH SIZE"
        if (cps /= huge(0)) then
          print *, "* Cells per side: ", cps
          write (*, '(1x, a, e10.3)') "* Domain size: ", domain_size
        end if
        print *, "* Global number of cells is ", mesh%topo%global_num_cells ! XXX: unknown at this point
        print *, "******************************************************************************"
        print *, "* RELAXATION FACTORS"
        write (*, '(1x, a, e10.3)') "* velocity: ", velocity_relax
        write (*, '(1x, a, e10.3)') "* pressure: ", pressure_relax
        print *, "******************************************************************************"
      end associate
      print *, "******************************************************************************"
      print *, "* REFERENCE VALUES"
      print *, "* Pressure      : ", run_options%reference_values%p_ref
      print *, "* Total Pressure: ", run_options%reference_values%p_total
      print *, "* Temperature   : ", run_options%reference_values%temp_ref
      print *, "* Density       : ", run_options%reference_values%dens_ref
      print *, "* Viscosity     : ", run_options%reference_values%visc_ref
      print *, "* Velocity      : ", run_options%reference_values%velo_ref
      print *, "* Length        : ", run_options%reference_values%len_ref
      print *, "* Reference cell: ", run_options%reference_values%pref_at_cell
      print *, "******************************************************************************"
    end if
    
  end subroutine
end submodule core_configuration
