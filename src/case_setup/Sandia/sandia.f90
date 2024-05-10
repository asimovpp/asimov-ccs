!> Program file for Sandia case
program sandia
#include "ccs_macros.inc"

  use petscvec
  use petscsys

  use core
  use ccs_base, only: mesh
  use case_config, only: case_name, &
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
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer, &
                      sync, is_root
  use parallel_types, only: parallel_environment
  use vec, only: create_vector, set_vector_location
  use petsctypes, only: vector_petsc
  use scalars, only: update_scalars
  use read_config, only: get_enable_cell_corrections, get_boundary_count, get_case_name, get_store_residuals
  use utils, only: set_size, initialise, update, exit_print, &
                   add_field_to_outputlist, get_field, add_field, &
                   set_is_field_solved, &
                   allocate_fluid_fields
  use boundary_conditions, only: read_bc_config, allocate_bc_arrays, set_bc_profile
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
  character(len=:), allocatable:: ccs_config_file  ! Config file for CCS
  character(len = ccs_string_len), dimension(:), allocatable:: variable_names  ! variable names for BC reading
  integer(ccs_int), dimension(:), allocatable:: variable_types              ! cell centred upwind, central, etc.

  type(ccs_options) :: run_options
  type(vector_spec) :: vec_properties

  type(field_spec):: field_properties
  class(field), pointer:: u, v, w, p, mf, viscosity, density
  class(field), pointer:: scalar_field

  integer(ccs_int):: n_boundaries

  integer(ccs_int):: irank  ! MPI rank ID
  integer(ccs_int):: isize  ! Size of MPI world

  integer(ccs_int):: timer_index_total
  integer(ccs_int):: timer_index_init
  integer(ccs_int):: i

  logical:: u_sol = .true.  ! Default equations to solve for LDC case
  logical:: v_sol = .true.
  logical:: w_sol = .true.
  logical:: p_sol = .true.

  logical:: store_residuals, enable_cell_corrections

  type(fluid):: flow_fields
  ! type(bc_profile), allocatable:: profile

  ! Launch MPI
  call initialise_parallel_environment(par_env)
  call timer_init()

  call get_config(par_env, run_options)
  call configure_parallelism(run_options, par_env, shared_env)
  irank = par_env%proc_id
  isize = par_env%num_procs

  call timer_register_start("Elapsed time", timer_index_total, is_total_time=.true.)

  call timer_register_start("Init time", timer_index_init)

  if (is_root(par_env)) print *, "Starting ", case_name, " case!"

  ! set solver and preconditioner info
  velocity_solver_method_name = "gmres"
  velocity_solver_precon_name = "bjacobi"
  pressure_solver_method_name = "cg"
  pressure_solver_precon_name = "gamg"

  ! Read mesh from .geo file
  call initialise_mesh(par_env, shared_env, run_options)

  ! Initialise fields
  if (irank == par_env%root) print *, "Initialise fields"

  ! Write gradients to solution file
  write_gradients = .false.

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

  ! Expect to find u, v, w, p, p_prime, scalar
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
  call get_field(flow_fields, "w", w)
  call get_field(flow_fields, "w", p)
  call get_field(flow_fields, "mf", mf)
  call get_field(flow_fields, "viscosity", viscosity)
  call get_field(flow_fields, "density", density)
  call get_field(flow_fields, "scalar", scalar_field)
  
  ! Add fields to output list
  call add_field_to_outputlist(u)
  call add_field_to_outputlist(v)
  call add_field_to_outputlist(w)
  call add_field_to_outputlist(p)
  call add_field_to_outputlist(scalar_field)

  ! Initialise velocity field
  if (irank == par_env%root) print *, "Initialise velocity field"
  call initialise_flow(flow_fields, get_init_flow, get_init_mass_flux)

  ! Solve using SIMPLE algorithm
  if (irank == par_env%root) print *, "Start SIMPLE"

  call activate_timestepping()
  call set_timestep(run_options%solve%dt)

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
  nullify(scalar_field)

  call timer_stop(timer_index_init)

  call run_solver(par_env, run_options, flow_fields)
  
  ! Clean-up

  call timer_stop(timer_index_total)

  call timer_print_all(par_env)
  call timer_export_csv(par_env)

  ! Finalise MPI
  call cleanup_parallel_environment(par_env)

contains

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
