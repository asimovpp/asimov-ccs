!v Program file for ScalarTransport case
!
!  @build mpi+petsc

program scalar_transport
#include "ccs_macros.inc"

  use petscvec
  use petscsys

  use core
  use case_config, only: case_name
  use constants, only: cell, face, ccsconfig, ccs_string_len, &
                       face_centred, cell_centred_central, cell_centred_upwind, &
                       ccs_split_type_low_high
  use kinds, only: ccs_real, ccs_int
  use types, only: field, field_spec, upwind_field, central_field, face_field, ccs_mesh, &
                   vector_spec, ccs_vector, field_ptr, fluid
  use fields, only: create_field, set_field_config_file, set_field_n_boundaries, set_field_name, &
                    set_field_type, set_field_vector_properties, set_field_store_residuals
  use parallel, only: initialise_parallel_environment, create_new_par_env, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments, &
                      is_root
  use parallel_types, only: parallel_environment
  use mesh_utils, only: build_mesh
  use meshing, only: get_global_num_cells, get_centre, count_neighbours, &
                     create_cell_locator, create_face_locator, create_neighbour_locator, &
                     get_local_index, get_boundary_status, get_face_normal, set_mesh_object, nullify_mesh_object, &
                     get_local_num_cells
  use vec, only: create_vector, set_vector_location, get_vector_data, restore_vector_data
  use scalars, only: update_scalars
  use utils, only: set_size, initialise, update, exit_print, add_field_to_outputlist, &
                   get_field, &
                   allocate_fluid_fields, dealloc_fluid_fields, &
                   get_scheme_name
  use boundary_conditions, only: read_bc_config, allocate_bc_arrays
  use read_config, only: get_store_residuals, get_variables, get_variable_types
  use timestepping, only: set_timestep, activate_timestepping, initialise_old_values, finalise_timestep

  implicit none

  class(parallel_environment), allocatable :: par_env
  class(parallel_environment), allocatable, target :: shared_env

  integer(ccs_int) :: irank ! MPI rank ID
  integer(ccs_int) :: isize ! Size of MPI world

  double precision :: start_time
  double precision :: init_time
  double precision :: end_time

  type(fluid) :: flow_fields

  real(ccs_real) :: L

  type(ccs_options) :: run_options
  
  ! Launch MPI
  call initialise_parallel_environment(par_env)

  call get_config(par_env, run_options)
  call configure_parallelism(run_options, par_env, shared_env)
  irank = par_env%proc_id
  isize = par_env%num_procs

  call timer(start_time)

  if (irank == par_env%root) print *, "Starting ", case_name, " case!"

  call initialise_mesh(par_env, shared_env, run_options)

  ! Initialise fields
  if (irank == par_env%root) print *, "Initialise fields"

  call initialise_fields(par_env, run_options, flow_fields)

  ! Initialise velocity field
  if (irank == par_env%root) print *, "Initialise flow field"
  call initialise_flow(par_env, run_options, flow_fields, get_init_flow, get_init_mass_flux)

  call activate_timestepping()
  call set_timestep(run_options%solve%dt)

  call timer(init_time)

  ! Solve using SIMPLE algorithm
  if (irank == par_env%root) print *, "Start scalar solver"

  ! Write out mesh and solution
  call run_solver(par_env, run_options, postproc_scalar, flow_fields)

  ! Clean-up
  call dealloc_fluid_fields(flow_fields)

  call timer(end_time)

  if (irank == par_env%root) then
    print *, "Init time: ", init_time - start_time
    print *, "Elapsed time: ", end_time - start_time
  end if

  ! Finalise MPI
  call nullify_mesh_object()
  call cleanup_parallel_environment(par_env)

contains

  pure subroutine get_init_flow(loc_p, field_name, init_val)

    use types, only: cell_locator

    type(cell_locator), intent(in) :: loc_p
    character(len=*), intent(in) :: field_name
    real(ccs_real), intent(inout) :: init_val

    real(ccs_real), dimension(3) :: c, r, x

    real(ccs_real) :: whisky

    c = L / 2

    call get_centre(loc_p, x)
    r = x - c

    if (any(r <= 0.0_ccs_real)) then
      whisky = 1.0_ccs_real
    else
      whisky = 0.0_ccs_real
    end if

    if (field_name == "whisky") then
      init_val = whisky
    else if (field_name == "water") then
      init_val = 1.0_ccs_real - whisky
    else
      init_val = init_val
    end if
    
  end subroutine get_init_flow

  pure subroutine get_init_mass_flux(loc_f, init_val)

    use constants, only: pi
    
    use types, only: face_locator

    type(face_locator), intent(in) :: loc_f
    real(ccs_real), intent(inout) :: init_val

    real(ccs_real), dimension(3) :: c, r, x
    real(ccs_real), dimension(3) :: face_normal, v
    real(ccs_real) :: rmag
    real(ccs_real) :: theta
    logical :: is_boundary
    
    c = L / 2

    call get_boundary_status(loc_f, is_boundary)
    if (is_boundary) then
      init_val = 0.0_ccs_real
    else
      call get_face_normal(loc_f, face_normal)
      call get_centre(loc_f, x)

      ! Create a rotating field that decays to zero on the boundaries
      r = x - c
      rmag = sqrt(sum(r(1:2)**2))
      theta = asin(abs(r(2)) / rmag)
      if (r(1) >= 0) then
        if (r(2) >= 0) then
          theta = theta + 0 * pi / 2
        else
          theta = 2 * pi - theta
        end if
      else
        if (r(2) >= 0) then
          theta = pi - theta
        else
          theta = theta + pi
        endif
      end if

      v(1) = -sin(theta)
      v(2) = cos(theta)
      v(3) = 0.0_ccs_real

      init_val = sum(v * face_normal)
    end if

  end subroutine get_init_mass_flux

  !> Subroutine to define case-specific postprocessing.
  subroutine postproc_scalar(par_env, flow_fields)

    class(parallel_environment), allocatable, intent(in) :: par_env
    type(fluid), intent(in) :: flow_fields

    ! All cases must define this, but if they don't require case-specific processing then simply
    ! make a no-op (use associate to silence unused variable compiler warnings)
    associate(foo => par_env, bar => flow_fields)
    end associate
    
  end subroutine postproc_scalar

end program scalar_transport
