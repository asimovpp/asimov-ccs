!> Program file for BackwardFacingStep case
program bfs
#include "ccs_macros.inc"

  use petscvec
  use petscsys

  use ccs_base, only: mesh
  use case_config, only: num_steps, num_iters, dt, domain_size, write_frequency, &
                         velocity_relax, pressure_relax, res_target, case_name, &
                         write_gradients, velocity_solver_method_name, velocity_solver_precon_name, &
                         pressure_solver_method_name, pressure_solver_precon_name, vertex_neighbours
  use constants, only: cell, face, ccsconfig, ccs_string_len, geoext, adiosconfig, ndim, &
                       field_u, field_v, field_w, field_p, field_p_prime, field_mf, field_viscosity, &
                       field_density, cell_centred_central, cell_centred_upwind, face_centred, &
                       ccs_split_type_shared, ccs_split_type_low_high, ccs_split_undefined
  use kinds, only: ccs_real, ccs_int, ccs_long
  use types, only: field, field_spec, upwind_field, central_field, face_field, ccs_mesh, &
                   vector_spec, ccs_vector, io_environment, io_process, &
                   field_ptr, fluid, fluid_solver_selector, bc_profile
  use fields, only: create_field, set_field_config_file, set_field_n_boundaries, set_field_name, &
                    set_field_type, set_field_vector_properties, set_field_store_residuals, set_field_enable_cell_corrections
  use fortran_yaml_c_interface, only: parse
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments, sync, &
                      create_new_par_env
  use parallel_types, only: parallel_environment
  use vec, only: create_vector, set_vector_location
  use petsctypes, only: vector_petsc
  use pv_coupling, only: solve_nonlinear
  use utils, only: set_size, initialise, update, exit_print, &
                   add_field_to_outputlist, get_field, add_field, &
                   get_fluid_solver_selector, set_fluid_solver_selector, &
                   allocate_fluid_fields
  use boundary_conditions, only: read_bc_config, allocate_bc_arrays, set_bc_profile
  use read_config, only: get_variables, get_boundary_count, get_case_name, get_store_residuals, get_enable_cell_corrections
  use timestepping, only: set_timestep, activate_timestepping, initialise_old_values
  use mesh_utils, only: read_mesh, write_mesh
  use meshing, only: set_mesh_object, nullify_mesh_object
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

  type(vector_spec) :: vec_properties

  type(field_spec) :: field_properties
  class(field), allocatable, target :: u, v, w, p, p_prime, mf, viscosity, density

  type(field_ptr), allocatable :: output_list(:)

  integer(ccs_int) :: n_boundaries

  integer(ccs_int) :: it_start, it_end
  integer(ccs_int) :: irank ! MPI rank ID
  integer(ccs_int) :: isize ! Size of MPI world

  double precision :: start_time
  double precision :: end_time

  logical :: u_sol = .true.  ! Default equations to solve for LDC case
  logical :: v_sol = .true.
  logical :: w_sol = .false.
  logical :: p_sol = .true.

  logical :: store_residuals, enable_cell_corrections

  integer(ccs_int) :: t          ! Timestep counter
  logical :: use_mpi_splitting

  type(fluid) :: flow_fields
  type(fluid_solver_selector) :: fluid_sol
  type(bc_profile), allocatable :: profile

  ! Launch MPI
  call initialise_parallel_environment(par_env)
  use_mpi_splitting = .true.
  call create_new_par_env(par_env, ccs_split_type_shared, use_mpi_splitting, shared_env)

  irank = par_env%proc_id
  isize = par_env%num_procs

  call read_command_line_arguments(par_env, case_name=case_name, in_dir=input_path)

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

  ! Hard coding to whether or not vertex neighbours are built
  vertex_neighbours = .false. ! set to .false. to avoid building

  ! Read mesh from .geo file
  if (irank == par_env%root) print *, "Reading mesh file"
  call read_mesh(par_env, shared_env, case_name, mesh)
  call set_mesh_object(mesh)

  ! Initialise fields
  if (irank == par_env%root) print *, "Initialise fields"

  ! Write gradients to solution file
  write_gradients = .true.

  ! Read boundary conditions
  if (irank == par_env%root) print *, "Read and allocate BCs"
  call get_boundary_count(ccs_config_file, n_boundaries)
  call get_store_residuals(ccs_config_file, store_residuals)
  call get_enable_cell_corrections(ccs_config_file, enable_cell_corrections)

  ! Create and initialise field vectors
  if (irank == par_env%root) print *, "Initialise field vectors"
  call initialise(vec_properties)

  call set_vector_location(cell, vec_properties)
  call set_size(par_env, mesh, vec_properties)

  call set_field_config_file(ccs_config_file, field_properties)
  call set_field_n_boundaries(n_boundaries, field_properties)
  call set_field_store_residuals(store_residuals, field_properties)
  call set_field_enable_cell_corrections(enable_cell_corrections, field_properties)

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
  call set_field_name("viscosity", field_properties)
  call create_field(field_properties, viscosity)
  call set_field_name("density", field_properties)
  call create_field(field_properties, density)

  call set_vector_location(face, vec_properties)
  call set_size(par_env, mesh, vec_properties)
  call set_field_vector_properties(vec_properties, field_properties)
  call set_field_type(face_centred, field_properties)
  call set_field_name("mf", field_properties)
  call create_field(field_properties, mf)

  ! Read and set BC profiles
  ! Read u componemt (1st column)
  call read_bc_profile(case_path // '.blasius.prf', 1, profile)
  profile%coordinates(:) = profile%coordinates(:) / mesh%geo%scalefactor
  profile%centre(:) = [ -4.0_ccs_real, 0.0_ccs_real, 0.5_ccs_real ] 

  ! Set to 3rd boundary condition (inlet)
  call set_bc_profile(u, profile, 3)

  ! Add fields to output list
  call add_field_to_outputlist(u, "u", output_list)
  call add_field_to_outputlist(v, "v", output_list)
  call add_field_to_outputlist(w, "w", output_list)
  call add_field_to_outputlist(p, "p", output_list)

  ! Initialise velocity field
  if (irank == par_env%root) print *, "Initialise velocity field"
  call initialise_flow(u, v, w, p, mf, viscosity, density)

  ! Solve using SIMPLE algorithm
  if (irank == par_env%root) print *, "Start SIMPLE"

  ! Write out mesh to file
  call write_mesh(par_env, case_path, mesh)

  ! Print the run configuration
  if (irank == par_env%root) then
    call print_configuration()
  end if

  call activate_timestepping()
  call set_timestep(dt)

  ! XXX: This should get incorporated as part of create_field subroutines
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

  do t = 1, num_steps
    call solve_nonlinear(par_env, mesh, it_start, it_end, res_target, &
                         fluid_sol, flow_fields)
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

  call nullify_mesh_object()
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
    print *, "* Global number of cells is ", mesh%topo%global_num_cells
    print *, "******************************************************************************"
    print *, "* RELAXATION FACTORS"
    write (*, '(1x,a,e10.3)') "* velocity: ", velocity_relax
    write (*, '(1x,a,e10.3)') "* pressure: ", pressure_relax
    print *, "******************************************************************************"

  end subroutine

  subroutine initialise_flow(u, v, w, p, mf, viscosity, density)

    use constants, only: insert_mode, ndim
    use types, only: vector_values, cell_locator, face_locator, neighbour_locator
    use meshing, only: create_cell_locator, get_global_index, count_neighbours, create_neighbour_locator, &
                       get_local_index, create_face_locator, get_local_index, get_face_normal, get_centre, &
                       get_local_num_cells
    use fv, only: calc_cell_coords
    use utils, only: clear_entries, set_mode, set_row, set_entry, set_values
    use vec, only: get_vector_data, restore_vector_data, create_vector_values

    ! Arguments
    class(field), intent(inout) :: u, v, w, p, mf, viscosity, density

    ! Local variables
    integer(ccs_int) :: n, count
    integer(ccs_int) :: n_local
    integer(ccs_int) :: index_p, global_index_p, index_f, index_nb
    real(ccs_real) :: u_val, v_val, w_val, p_val
    type(cell_locator) :: loc_p
    type(face_locator) :: loc_f
    type(neighbour_locator) :: loc_nb
    type(vector_values) :: u_vals, v_vals, w_vals, p_vals
    real(ccs_real), dimension(:), pointer :: mf_data, viscosity_data, density_data

    real(ccs_real), dimension(ndim) :: x_p, x_f
    real(ccs_real), dimension(ndim) :: face_normal

    integer(ccs_int) :: nnb
    integer(ccs_int) :: j

    ! Set alias
    call get_local_num_cells(n_local)

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
      call create_cell_locator(index_p, loc_p)
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
          mf_data(index_f) = 0.0_ccs_real
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
    call update(w%values)
    call update(p%values)
    call update(mf%values)
    call update(viscosity%values)
    call update(density%values)

  end subroutine initialise_flow

  subroutine read_bc_profile(filename, variable_id, profile)
    
    character(len=*), intent(in) :: filename
    integer(ccs_int), intent(in) :: variable_id
    type(bc_profile), allocatable, intent(out) :: profile

    real(ccs_real), allocatable, dimension(:) :: tmp_values
    real(ccs_real) :: tmp_coord
    character(len=128) :: header_string, tmp
    integer(ccs_int) :: num_field, i
    integer :: io_err, unit_io

    allocate(profile)

    allocate(profile%centre(3))
    allocate(profile%values(0))
    allocate(profile%coordinates(0))

    open(newunit=unit_io, file=trim(filename), status='old', action='read')

    read(unit_io, *)                      ! ignore profile type
    read(unit_io, *) tmp, profile%centre  ! read centre
    read(unit_io, *)                      ! ignore tolerance
    read(unit_io, *)                      ! ignore scaling
    read(unit_io, '(A)') header_string

    ! Count the number of fields in file
    num_field = -1
    do i=1, len(header_string)
      if (header_string(i:i) == ',') then
        num_field = num_field + 1
      end if
    end do

    allocate(tmp_values(num_field))

    ! Read file profile table
    do while (.true.)

      read(unit_io, *, iostat=io_err) tmp_coord, tmp_values
      if (io_err /= 0) then
        exit
      end if

      profile%values = [ profile%values, tmp_values(variable_id) ]
      profile%coordinates = [ profile%coordinates, tmp_coord ]
    end do

  end subroutine read_bc_profile

end program bfs
