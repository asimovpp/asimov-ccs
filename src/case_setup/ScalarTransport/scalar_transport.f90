!v Program file for ScalarTransport case
!
!  @build mpi+petsc

program scalar_transport
#include "ccs_macros.inc"

  use petscvec
  use petscsys

  use case_config, only: num_steps, num_iters, cps, domain_size, case_name, &
                         velocity_relax, pressure_relax, res_target, &
                         write_gradients, velocity_solver_method_name, velocity_solver_precon_name, &
                         pressure_solver_method_name, pressure_solver_precon_name, &
                         dt, write_frequency
  use constants, only: cell, face, ccsconfig, ccs_string_len, field_u, field_v, &
                       field_w, field_p, field_p_prime, field_mf, &
                       cell_centred_central, cell_centred_upwind, face_centred
  use kinds, only: ccs_real, ccs_int
  use types, only: field, field_spec, upwind_field, central_field, face_field, ccs_mesh, &
                   vector_spec, ccs_vector, field_ptr, fluid, fluid_solver_selector
  use fields, only: create_field, set_field_config_file, set_field_n_boundaries, set_field_name, &
                    set_field_type, set_field_vector_properties, set_field_store_residuals
  use fortran_yaml_c_interface, only: parse
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments, sync
  use parallel_types, only: parallel_environment
  use mesh_utils, only: build_mesh, write_mesh, build_square_mesh
  use meshing, only: get_global_num_cells, get_centre, count_neighbours, &
                     create_cell_locator, create_face_locator, create_neighbour_locator, &
                     get_local_index, get_boundary_status, get_face_normal
  use vec, only: create_vector, set_vector_location
  use petsctypes, only: vector_petsc
  use scalars, only: calculate_scalars
  use utils, only: set_size, initialise, update, exit_print, add_field_to_outputlist, &
                   get_field, set_field, get_fluid_solver_selector, set_fluid_solver_selector, &
                   allocate_fluid_fields, dealloc_fluid_fields
  use boundary_conditions, only: read_bc_config, allocate_bc_arrays
  use read_config, only: get_bc_variables, get_boundary_count, get_store_residuals
  use io_visualisation, only: write_solution
  use timestepping, only: set_timestep, activate_timestepping, initialise_old_values

  implicit none

  class(parallel_environment), allocatable :: par_env
  character(len=:), allocatable :: input_path  ! Path to input directory
  character(len=:), allocatable :: case_path  ! Path to input directory with case name appended
  character(len=:), allocatable :: ccs_config_file ! Config file for CCS
  character(len=ccs_string_len), dimension(:), allocatable :: variable_names  ! variable names for BC reading

  type(ccs_mesh) :: mesh
  type(vector_spec) :: vec_properties

  type(field_spec) :: field_properties
  class(field), allocatable, target :: whisky, water, mf

  type(field_ptr), allocatable :: output_list(:)

  integer(ccs_int) :: n_boundaries

  integer(ccs_int) :: it_start, it_end
  integer(ccs_int) :: irank ! MPI rank ID
  integer(ccs_int) :: isize ! Size of MPI world

  double precision :: start_time
  double precision :: init_time
  double precision :: end_time

  logical :: u_sol = .true.  ! Default equations to solve for LDC case
  logical :: v_sol = .true.
  logical :: w_sol = .true.
  logical :: p_sol = .true.

  logical :: store_residuals

  type(fluid) :: flow_fields
  type(fluid_solver_selector) :: fluid_sol

  real(ccs_real) :: L
  integer(ccs_int) :: t

  integer(ccs_int) :: scalar_index ! ID for defining scalars
  
  ! Launch MPI
  call initialise_parallel_environment(par_env)

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
  L = 1.0_ccs_real
  mesh = build_mesh(par_env, cps, cps, cps, L)   ! 3-D mesh

  ! Initialise fields
  if (irank == par_env%root) print *, "Initialise fields"

  ! Write gradients to solution file
  write_gradients = .false.

  ! Create and initialise field vectors
  call get_boundary_count(ccs_config_file, n_boundaries)
  call get_bc_variables(ccs_config_file, variable_names)
  call get_store_residuals(ccs_config_file, store_residuals)

  call initialise(vec_properties)

  call set_vector_location(cell, vec_properties)
  call set_size(par_env, mesh, vec_properties)

  call set_field_config_file(ccs_config_file, field_properties)
  call set_field_n_boundaries(n_boundaries, field_properties)
  call set_field_store_residuals(store_residuals, field_properties)

  call set_field_vector_properties(vec_properties, field_properties)
  call set_field_type(cell_centred_upwind, field_properties)
  call set_field_name("whisky", field_properties)
  call create_field(field_properties, whisky)
  call set_field_name("water", field_properties)
  call create_field(field_properties, water)

  call set_vector_location(face, vec_properties)
  call set_size(par_env, mesh, vec_properties)
  call set_field_vector_properties(vec_properties, field_properties)
  call set_field_type(face_centred, field_properties)
  call set_field_name("mf", field_properties)
  call create_field(field_properties, mf)

  ! Add fields to output list
  allocate (output_list(2))
  call add_field_to_outputlist(whisky, "whisky", output_list)
  call add_field_to_outputlist(water, "water", output_list)

  ! Initialise velocity field
  if (irank == par_env%root) print *, "Initialise flow field"
  call initialise_case(mesh, whisky, water, mf)
  call update(whisky%values)
  call update(water%values)
  call update(mf%values)

  ! ! XXX: This should get incorporated as part of create_field subroutines
  ! call set_fluid_solver_selector(field_u, u_sol, fluid_sol)
  ! call set_fluid_solver_selector(field_v, v_sol, fluid_sol)
  ! call set_fluid_solver_selector(field_w, w_sol, fluid_sol)
  ! call set_fluid_solver_selector(field_p, p_sol, fluid_sol)
  call allocate_fluid_fields(3, flow_fields)
  call set_field(1, field_mf, mf, flow_fields)
  scalar_index = maxval((/ field_u, field_v, field_w, field_p, field_p_prime, field_mf /)) + 1
  call set_field(2, scalar_index, whisky, flow_fields)
  scalar_index = scalar_index + 1
  call set_field(3, scalar_index, water, flow_fields)

  if (irank == par_env%root) then
    call print_configuration()
  end if

  call activate_timestepping()
  call set_timestep(dt)

  call timer(init_time)
  ! Solve using SIMPLE algorithm
  if (irank == par_env%root) print *, "Start scalar solver"

  ! ! Write out mesh and solution
  call write_mesh(par_env, case_path, mesh)
  call write_solution(par_env, case_path, mesh, output_list, 0, num_steps, dt)
  do t = 1, num_steps
    call calculate_scalars(par_env, mesh, flow_fields)
    if (par_env%proc_id == par_env%root) then
      print *, "TIME = ", t, " / ", num_steps
    end if

    if ((t == 1) .or. (t == num_steps) .or. (mod(t, write_frequency) == 0)) then
      call write_solution(par_env, case_path, mesh, output_list, t, num_steps, dt)
    end if
  end do

  ! Clean-up
  call dealloc_fluid_fields(flow_fields)
  deallocate (whisky)
  deallocate (water)
  deallocate (mf)
  deallocate (output_list)

  call timer(end_time)

  if (irank == par_env%root) then
    print *, "Init time: ", init_time - start_time
    print *, "Elapsed time: ", end_time - start_time
  end if

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

    integer(ccs_int) :: global_num_cells

    call get_global_num_cells(mesh, global_num_cells)

    ! XXX: this should eventually be replaced by something nicely formatted that uses "write"
    print *, " "
    print *, "******************************************************************************"
    print *, "* Solving the ", case_name, " case"
    print *, "******************************************************************************"
    print *, " "
    print *, "******************************************************************************"
    print *, "* SIMULATION LENGTH"
    print *, "* Running for ", num_iters, "iterations"
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

  subroutine initialise_case(mesh, whisky, water, mf)

    use constants, only: add_mode
    use types, only: vector_values, cell_locator, neighbour_locator, face_locator
    use meshing, only: create_cell_locator, get_global_index, get_local_num_cells
    use fv, only: calc_cell_coords
    use utils, only: clear_entries, set_mode, set_row, set_entry, set_values
    use vec, only: get_vector_data, restore_vector_data, create_vector_values

    ! Arguments
    class(ccs_mesh), intent(in) :: mesh
    class(field), intent(inout) :: whisky, water, mf

    ! Local variables
    integer(ccs_int) :: row, col
    integer(ccs_int) :: index_p, global_index_p, n_local
    real(ccs_real) :: whisky_val, water_val
    type(cell_locator) :: loc_p
    type(face_locator) :: loc_f
    type(neighbour_locator) :: loc_nb
    type(vector_values) :: whisky_vals, water_vals
    real(ccs_real), dimension(:), pointer :: mf_data

    integer(ccs_int) :: j, nnb
    integer(ccs_int) :: index_nb, index_f
    logical :: is_boundary
    
    real(ccs_real), dimension(3) :: x, r, c, face_normal
    real(ccs_real), dimension(3) :: v
    real(ccs_real) :: theta, rmag
    
    ! Set alias
    call get_local_num_cells(mesh, n_local)

    call create_vector_values(n_local, whisky_vals)
    call create_vector_values(n_local, water_vals)
    call set_mode(add_mode, whisky_vals)
    call set_mode(add_mode, water_vals)

    call get_vector_data(mf%values, mf_data)
    
    ! Coordinates of domain centre
    c = L / 2
    
    ! Set initial values for velocity fields
    do index_p = 1, n_local
      call create_cell_locator(mesh, index_p, loc_p)
      call get_global_index(loc_p, global_index_p)

      ! Get domain centre offset
      call get_centre(loc_p, x)
      r = c - x

      if (any(r <= 0.0_ccs_real)) then
        whisky_val = 1.0_ccs_real
      else
        whisky_val = 0.0_ccs_real
      end if
      water_val = 1.0_ccs_real - whisky_val

      call set_row(global_index_p, whisky_vals)
      call set_entry(whisky_val, whisky_vals)
      call set_row(global_index_p, water_vals)
      call set_entry(water_val, water_vals)

      call count_neighbours(loc_p, nnb)
      do j = 1, nnb
         call create_face_locator(mesh, index_p, j, loc_f)
         call get_local_index(loc_f, index_f)
         call get_boundary_status(loc_f, is_boundary)

         if (.not. is_boundary) then
            call create_neighbour_locator(loc_p, j, loc_nb)
            call get_local_index(loc_nb, index_nb)

            if (index_nb > index_p) then
               call create_face_locator(mesh, index_p, j, loc_f)
               call get_face_normal(loc_f, face_normal)
               call get_centre(loc_f, x)

               ! Create a rotating field that decays to zero on the boundaries
               r = c - x
               rmag = sqrt(sum(r(1:2)**2))
               theta = asin(r(1) / rmag)

               v(1) = -sin(theta)
               v(2) = cos(theta)
               v(3) = 0.0_ccs_real
               
               mf_data(index_f) = sum(v * face_normal)
               mf_data(index_f) = 0.0 
            end if
         else
            mf_data(index_f) = 0.0_ccs_real
         end if
      end do
      
    end do

    call set_values(whisky_vals, whisky%values)
    call set_values(water_vals, water%values)

    deallocate (whisky_vals%global_indices)
    deallocate (water_vals%global_indices)
    deallocate (whisky_vals%values)
    deallocate (water_vals%values)

    call restore_vector_data(mf%values, mf_data)

  end subroutine initialise_case

end program scalar_transport
