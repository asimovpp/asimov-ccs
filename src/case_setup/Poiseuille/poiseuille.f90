!> Program file for Poiseuille case
program poiseuille
#include "ccs_macros.inc"

  use petscvec
  use petscsys

  use case_config, only: num_steps, num_iters, dt, cps, domain_size, write_frequency, &
                         velocity_relax, pressure_relax, res_target, case_name, &
                         write_gradients, velocity_solver_method_name, velocity_solver_precon_name, &
                         pressure_solver_method_name, pressure_solver_precon_name
  use constants, only: ccs_split_type_shared, ccs_split_type_low_high, ccs_split_undefined
  use constants, only: cell, face, ccsconfig, ccs_string_len, geoext, adiosconfig, ndim, &
                       field_u, field_v, field_w, field_p, field_p_prime, field_mf, &
                       cell_centred_central, cell_centred_upwind, face_centred
  use kinds, only: ccs_real, ccs_int, ccs_long
  use types, only: field, field_spec, upwind_field, central_field, face_field, ccs_mesh, &
                   vector_spec, ccs_vector, io_environment, io_process, &
                   field_ptr, fluid, fluid_solver_selector
  use fields, only: create_field, set_field_config_file, set_field_n_boundaries, set_field_name, &
       set_field_type, set_field_vector_properties
  use fortran_yaml_c_interface, only: parse
  use parallel, only: initialise_parallel_environment, create_new_par_env, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments, sync
  use parallel_types, only: parallel_environment
  use vec, only: create_vector, set_vector_location
  use petsctypes, only: vector_petsc
  use pv_coupling, only: solve_nonlinear
  use utils, only: set_size, initialise, update, exit_print, &
                   calc_kinetic_energy, calc_enstrophy, &
                   add_field_to_outputlist, get_field, set_field, &
                   get_fluid_solver_selector, set_fluid_solver_selector, &
                   allocate_fluid_fields
  use boundary_conditions, only: read_bc_config, allocate_bc_arrays
  use read_config, only: get_variables, get_boundary_count, get_case_name
  use timestepping, only: set_timestep, activate_timestepping, initialise_old_values
  use mesh_utils, only: read_mesh, build_square_mesh, write_mesh
  use partitioning, only: compute_partitioner_input, &
                          partition_kway, compute_connectivity
  use io_visualisation, only: write_solution
  use fv, only: update_gradient
  use utils, only: str

  implicit none

  class(parallel_environment), allocatable :: par_env
  class(parallel_environment), allocatable :: shared_env
  character(len=:), allocatable :: input_path  ! Path to input directory
  character(len=:), allocatable :: case_path  ! Path to input directory with case name appended
  character(len=:), allocatable :: ccs_config_file ! Config file for CCS
  character(len=ccs_string_len), dimension(:), allocatable :: variable_names  ! variable names for BC reading

  type(ccs_mesh) :: mesh

  type(vector_spec) :: vec_properties

  type(field_spec) :: field_properties
  class(field), allocatable, target :: u, v, w, p, p_prime, mf

  type(field_ptr), allocatable :: output_list(:)

  integer(ccs_int) :: n_boundaries

  integer(ccs_int) :: it_start, it_end
  integer(ccs_int) :: irank ! MPI rank ID
  integer(ccs_int) :: isize ! Size of MPI world

  double precision :: start_time
  double precision :: end_time

  logical :: u_sol = .true.  ! Default equations to solve for LDC case
  logical :: v_sol = .true.
  logical :: w_sol = .true.
  logical :: p_sol = .true.

  integer(ccs_int) :: t          ! Timestep counter
  logical :: use_mpi_splitting

  type(fluid) :: flow_fields
  type(fluid_solver_selector) :: fluid_sol

  ! Launch MPI
  call initialise_parallel_environment(par_env)
  use_mpi_splitting = .false.
  call create_new_par_env(par_env, ccs_split_type_low_high, use_mpi_splitting, shared_env)

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

  ! Set start and end iteration numbers (read from input file)
  it_start = 1
  it_end = num_iters

  ! If cps is no longer the default value, it has been set explicity and
  ! the mesh generator is invoked...
  if (cps /= huge(0)) then
    ! Create a cubic mesh
    if (irank == par_env%root) print *, "Building mesh"
    mesh = build_square_mesh(par_env, shared_env, cps, domain_size)
  else
    if (irank == par_env%root) print *, "Reading mesh file"
    call read_mesh(par_env, shared_env, case_name, mesh)
  end if

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

  call set_field_config_file(ccs_config_file, field_properties)
  call set_field_n_boundaries(n_boundaries, field_properties)

  call set_field_vector_properties(vec_properties, field_properties)
  call set_field_type(cell_centred_upwind, field_properties)
  call set_field_name("u", field_properties)
  call create_field(field_properties, u)
  call set_field_name("v", field_properties)
  call create_field(field_properties, v)
  call set_field_name("w", field_properties)
  call create_field(field_properties, w)

  call set_field_type(cell_centred_central, field_properties)
  call set_field_name("p", field_properties)
  call create_field(field_properties, p)
  call set_field_name("p_prime", field_properties)
  call create_field(field_properties, p_prime)

  call set_vector_location(face, vec_properties)
  call set_size(par_env, mesh, vec_properties)
  call set_field_vector_properties(vec_properties, field_properties)
  call set_field_type(face_centred, field_properties)
  call set_field_name("mf", field_properties)
  call create_field(field_properties, mf)

  ! Add fields to output list
  call add_field_to_outputlist(u, "u", output_list)
  call add_field_to_outputlist(v, "v", output_list)
  call add_field_to_outputlist(w, "w", output_list)
  call add_field_to_outputlist(p, "p", output_list)

  call activate_timestepping()
  call set_timestep(dt)

  ! Initialise velocity field
  if (irank == par_env%root) print *, "Initialise velocity field"
  call initialise_flow(mesh, u, v, w, p, mf)
  call calc_kinetic_energy(par_env, mesh, u, v, w)
  call calc_enstrophy(par_env, mesh, u, v, w)

  ! Solve using SIMPLE algorithm
  if (irank == par_env%root) print *, "Start SIMPLE"
  call calc_kinetic_energy(par_env, mesh, u, v, w)
  call calc_enstrophy(par_env, mesh, u, v, w)

  ! Write out mesh to file
  call write_mesh(par_env, case_path, mesh)

  ! Print the run configuration
  if (irank == par_env%root) then
    call print_configuration()
  end if

  ! XXX: This should get incorporated as part of create_field subroutines
  call set_fluid_solver_selector(field_u, u_sol, fluid_sol)
  call set_fluid_solver_selector(field_v, v_sol, fluid_sol)
  call set_fluid_solver_selector(field_w, w_sol, fluid_sol)
  call set_fluid_solver_selector(field_p, p_sol, fluid_sol)
  call allocate_fluid_fields(6, flow_fields)
  call set_field(1, field_u, u, flow_fields)
  call set_field(2, field_v, v, flow_fields)
  call set_field(3, field_w, w, flow_fields)
  call set_field(4, field_p, p, flow_fields)
  call set_field(5, field_p_prime, p_prime, flow_fields)
  call set_field(6, field_mf, mf, flow_fields)

  do t = 1, num_steps
    call solve_nonlinear(par_env, mesh, it_start, it_end, res_target, &
                         fluid_sol, flow_fields)
    call calc_kinetic_energy(par_env, mesh, u, v, w)
    call update_gradient(mesh, u)
    call update_gradient(mesh, v)
    call update_gradient(mesh, w)
    call calc_enstrophy(par_env, mesh, u, v, w)
    if (par_env%proc_id == par_env%root) then
      print *, "TIME = ", t
    end if

    if ((t == 1) .or. (t == num_steps) .or. (mod(t, write_frequency) == 0)) then
      call write_solution(par_env, case_path, mesh, output_list, t, num_steps, dt)
    end if
  end do

  ! Clean-up
  deallocate (u)
  deallocate (v)
  deallocate (w)
  deallocate (p)
  deallocate (p_prime)
  deallocate (output_list)

  call timer(end_time)

  if (irank == par_env%root) then
    print *, "Elapsed time: ", end_time - start_time
  end if

  ! Finalise MPI
  call cleanup_parallel_environment(par_env)

contains

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

    if (cps == huge(0)) then ! cps was not set on the command line
      call get_value(config_file, 'cps', cps)
      if (cps == huge(0)) then
        call error_abort("No value assigned to cps.")
      end if
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

    call get_relaxation_factors(config_file, u_relax=velocity_relax, p_relax=pressure_relax)
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
    write (*, '(1x,a,e10.3)') "* Time step size: ", dt
    print *, "******************************************************************************"
    print *, "* MESH SIZE"
    if (cps /= huge(0)) then
      print *, "* Cells per side: ", cps
      write (*, '(1x,a,e10.3)') "* Domain size: ", domain_size
    end if
    print *, "* Global number of cells is ", mesh%topo%global_num_cells
    print *, "******************************************************************************"
    print *, "* RELAXATION FACTORS"
    write (*, '(1x,a,e10.3)') "* velocity: ", velocity_relax
    write (*, '(1x,a,e10.3)') "* pressure: ", pressure_relax
    print *, "******************************************************************************"

  end subroutine

  subroutine initialise_flow(mesh, u, v, w, p, mf)

    use constants, only: insert_mode, ndim
    use types, only: vector_values, cell_locator, face_locator, neighbour_locator
    use meshing, only: create_cell_locator, get_global_index, count_neighbours, create_neighbour_locator, &
                       get_local_index, create_face_locator, get_local_index, get_face_normal, get_centre, &
                       get_local_num_cells
    use fv, only: calc_cell_coords
    use utils, only: clear_entries, set_mode, set_row, set_entry, set_values
    use vec, only: get_vector_data, restore_vector_data, create_vector_values

    ! Arguments
    class(ccs_mesh), intent(in) :: mesh
    class(field), intent(inout) :: u, v, w, p, mf

    ! Local variables
    integer(ccs_int) :: n, count
    integer(ccs_int) :: n_local
    integer(ccs_int) :: index_p, global_index_p, index_f, index_nb
    real(ccs_real) :: u_val, v_val, w_val, p_val
    type(cell_locator) :: loc_p
    type(face_locator) :: loc_f
    type(neighbour_locator) :: loc_nb
    type(vector_values) :: u_vals, v_vals, w_vals, p_vals
    real(ccs_real), dimension(:), pointer :: mf_data

    real(ccs_real), dimension(ndim) :: x_p, x_f
    real(ccs_real), dimension(ndim) :: face_normal

    integer(ccs_int) :: nnb
    integer(ccs_int) :: j

    ! Set alias
    call get_local_num_cells(mesh, n_local)

    call create_vector_values(n_local, u_vals)
    call create_vector_values(n_local, v_vals)
    call create_vector_values(n_local, w_vals)
    call create_vector_values(n_local, p_vals)
    call set_mode(insert_mode, u_vals)
    call set_mode(insert_mode, v_vals)
    call set_mode(insert_mode, w_vals)
    call set_mode(insert_mode, p_vals)

    ! Set initial values for velocity fields
    do index_p = 1, n_local
      call create_cell_locator(mesh, index_p, loc_p)
      call get_global_index(loc_p, global_index_p)

      call get_centre(loc_p, x_p)

      u_val = 0.0_ccs_real 
      v_val = 0.0_ccs_real 
      w_val = 0.0_ccs_real
      p_val = 0.0_ccs_real

      call set_row(global_index_p, u_vals)
      call set_entry(u_val, u_vals)
      call set_row(global_index_p, v_vals)
      call set_entry(v_val, v_vals)
      call set_row(global_index_p, w_vals)
      call set_entry(w_val, w_vals)
      call set_row(global_index_p, p_vals)
      call set_entry(p_val, p_vals)
    end do

    call set_values(u_vals, u%values)
    call set_values(v_vals, v%values)
    call set_values(w_vals, w%values)
    call set_values(p_vals, p%values)

    deallocate (u_vals%global_indices)
    deallocate (v_vals%global_indices)
    deallocate (w_vals%global_indices)
    deallocate (p_vals%global_indices)
    deallocate (u_vals%values)
    deallocate (v_vals%values)
    deallocate (w_vals%values)
    deallocate (p_vals%values)

    call get_vector_data(mf%values, mf_data)

    count = 0
    n = 0

    ! Loop over local cells and faces
    call get_local_num_cells(mesh, n_local)
    do index_p = 1, n_local

      call create_cell_locator(mesh, index_p, loc_p)
      call count_neighbours(loc_p, nnb)
      do j = 1, nnb

        call create_neighbour_locator(loc_p, j, loc_nb)
        call get_local_index(loc_nb, index_nb)

        ! if neighbour index is greater than previous face index
        if (index_nb > index_p) then ! XXX: abstract this test

          call create_face_locator(mesh, index_p, j, loc_f)
          call get_local_index(loc_f, index_f)
          call get_face_normal(loc_f, face_normal)
          call get_centre(loc_f, x_f)

          ! compute initial value based on current face coordinates
          mf_data(index_f) = 0.0_ccs_real * face_normal(1)
        end if

      end do
    end do

    call restore_vector_data(mf%values, mf_data)

    call update(u%values)
    call update(v%values)
    call update(w%values)
    call update(p%values)
    call update(mf%values)

  end subroutine initialise_flow

end program poiseuille
