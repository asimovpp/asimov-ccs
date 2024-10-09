!> Program file for TaylorGreenVortex case
program tgv
#include "ccs_macros.inc"

  use petscsys
  use petscvec

  use core
  use boundary_conditions, only: read_bc_config, allocate_bc_arrays
  use case_config, only: case_name, write_gradients
  use constants, only: cell, face, ccsconfig, ccs_string_len, geoext, adiosconfig, ndim, &
                       cell_centred_central, cell_centred_upwind, face_centred
  use constants, only: ccs_split_type_shared, ccs_split_type_low_high, ccs_split_undefined
  use meshing, only: set_mesh_object, nullify_mesh_object
  use fields, only: create_field, set_field_config_file, set_field_n_boundaries, set_field_name, &
                    set_field_type, set_field_vector_properties, set_field_store_residuals, set_field_enable_cell_corrections
  use fortran_yaml_c_interface, only: parse
  use fv, only: update_gradient
  use kinds, only: ccs_real, ccs_int, ccs_long
  use mesh_utils, only: read_mesh, build_mesh
  use parallel, only: initialise_parallel_environment, &
                      create_new_par_env, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments, sync, is_root
  use parallel_types, only: parallel_environment
  use partitioning, only: compute_partitioner_input, &
                          partition_kway, compute_connectivity
  use petsctypes, only: vector_petsc
  use read_config, only: get_variables, get_case_name, get_store_residuals, get_enable_cell_corrections, get_variable_types
  use timestepping, only: set_timestep, activate_timestepping, initialise_old_values
  use types, only: field, field_spec, upwind_field, central_field, face_field, ccs_mesh, &
                   vector_spec, ccs_vector, io_environment, io_process, &
                   field_ptr, fluid
  use utils, only: set_size, initialise, update, exit_print, &
                   calc_kinetic_energy, calc_enstrophy, &
                   add_field_to_outputlist, get_field, add_field, &
                   set_is_field_solved, &
                   allocate_fluid_fields, str, debug_print, &
                   dealloc_fluid_fields
  use vec, only: create_vector, set_vector_location
  use timers, only: timer_init, timer_register_start, timer_register, timer_start, timer_stop, timer_print, &
                    timer_get_time, timer_print_all, timer_export_csv

  implicit none

  class(parallel_environment), allocatable:: par_env
  class(parallel_environment), allocatable:: shared_env

  type(ccs_options) :: run_options

  integer(ccs_int):: irank  ! MPI rank ID
  integer(ccs_int):: isize  ! Size of MPI world

  integer(ccs_int):: timer_index_total
  integer(ccs_int):: timer_index_init
  integer(ccs_int):: timer_index_io_sol
  integer(ccs_int):: timer_index_sol

  double precision:: sol_time, io_time

  type(fluid):: flow_fields

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

  call initialise_mesh(par_env, shared_env, run_options)

  ! Initialise fields
  if (irank == par_env%root) print *, "Initialise fields"

  ! Write gradients to solution file
  write_gradients = .true.

  call initialise_fields(par_env, run_options, flow_fields)
  
  call activate_timestepping()
  call set_timestep(run_options%solve%dt)

  ! Initialise velocity field
  if (irank == par_env%root) print *, "Initialise velocity field"
  call initialise_flow(par_env, run_options, flow_fields, get_init_flow, get_init_mass_flux)

  call timer_stop(timer_index_init)
  call timer_register("I/O time for solution", timer_index_io_sol)
  call timer_register("Solver time inc I/O", timer_index_sol)

  call run_solver(par_env, run_options, eval_sources, postproc_tgv, flow_fields)

  call timer_stop(timer_index_total)

  call timer_print_all(par_env)
  call timer_export_csv(par_env)

  call timer_get_time(timer_index_sol, sol_time)
  call timer_get_time(timer_index_io_sol, io_time)
  if (irank == par_env%root) then
    write(*,'(A30, F10.4, A)') "Solver time no I/O:", sol_time - io_time, " s"
    write(*,'(A30, F10.4, A)') "Average time/step (no I/O):", (sol_time - io_time) / run_options%solve%num_steps, " s"
  end if

  call dealloc_fluid_fields(flow_fields)
  call nullify_mesh_object()
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

    select case (field_name)
    case ("u")
      call get_centre(loc_p, x_p)
      init_val = sin(x_p(1)) * cos(x_p(2)) * cos(x_p(3))
    case ("v")
      call get_centre(loc_p, x_p)
      init_val = -cos(x_p(1)) * sin(x_p(2)) * cos(x_p(3))
    case default
      init_val = init_val + 0.0_ccs_real ! Keep the default value
    end select
    
  end subroutine get_init_flow
  
  pure subroutine get_init_mass_flux(loc_f, init_val)
    use types, only: face_locator
    use meshing, only: get_face_normal, get_centre
    type(face_locator), intent(in) :: loc_f
    real(ccs_real), intent(inout) :: init_val
    real(ccs_real), dimension(ndim) :: x_f
    real(ccs_real), dimension(ndim) :: face_normal

    call get_face_normal(loc_f, face_normal)
    call get_centre(loc_f, x_f)
  
    init_val = sin(x_f(1)) * cos(x_f(2)) * cos(x_f(3)) * face_normal(1) &
             - cos(x_f(1)) * sin(x_f(2)) * cos(x_f(3)) * face_normal(2)

  end subroutine get_init_mass_flux

  subroutine postproc_tgv(par_env, flow_fields)

    class(parallel_environment), allocatable, intent(in) :: par_env
    type(fluid), intent(in) :: flow_fields

    class(field), pointer:: u, v, w

    call get_field(flow_fields, "u", u)
    call get_field(flow_fields, "v", v)
    call get_field(flow_fields, "w", w)
    call calc_kinetic_energy(par_env, u, v, w)
    call calc_enstrophy(par_env, u, v, w)
    nullify(u)
    nullify(v)
    nullify(w)

  end subroutine postproc_tgv

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

end program tgv
