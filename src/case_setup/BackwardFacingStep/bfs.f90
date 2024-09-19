!> Program file for BackwardFacingStep case
program bfs
#include "ccs_macros.inc"

  use petscvec
  use petscsys

  use core
  use ccs_base, only: mesh
  use case_config, only: case_name, write_gradients
  use constants, only: cell, face, ccsconfig, ccs_string_len, geoext, adiosconfig, ndim, &
                       cell_centred_central, cell_centred_upwind, face_centred, &
                       ccs_split_type_shared, ccs_split_type_low_high, ccs_split_undefined
  use kinds, only: ccs_real, ccs_int, ccs_long
  use types, only: field, field_spec, upwind_field, central_field, face_field, ccs_mesh, &
                   vector_spec, ccs_vector, io_environment, io_process, &
                   field_ptr, fluid, bc_profile
  use fields, only: create_field, set_field_config_file, set_field_n_boundaries, set_field_name, &
                    set_field_type, set_field_vector_properties, set_field_store_residuals, set_field_enable_cell_corrections
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments, sync, &
                      create_new_par_env, is_root
  use parallel_types, only: parallel_environment
  use vec, only: create_vector, set_vector_location
  use utils, only: set_size, initialise, update, exit_print, &
                   add_field_to_outputlist, get_field, add_field, &
                   set_is_field_solved, &
                   allocate_fluid_fields
  use boundary_conditions, only: read_bc_config, allocate_bc_arrays, set_bc_profile
  use read_config, only: get_variables, get_case_name, &
                         get_store_residuals, get_enable_cell_corrections, get_variable_types
  use timestepping, only: set_timestep, activate_timestepping, initialise_old_values
  use mesh_utils, only: read_mesh
  use meshing, only: set_mesh_object, nullify_mesh_object
  use partitioning, only: compute_partitioner_input, &
                          partition_kway, compute_connectivity
  use fv, only: update_gradient
  use utils, only: str
  use timers, only: timer_print_all, timer_export_csv

  implicit none

  class(parallel_environment), allocatable :: par_env
  class(parallel_environment), allocatable :: shared_env
  character(len=:), allocatable :: case_path  ! Path to input directory with case name appended

  class(field), pointer :: u

  integer(ccs_int) :: irank ! MPI rank ID
  integer(ccs_int) :: isize ! Size of MPI world

  double precision :: start_time
  double precision :: end_time

  type(fluid) :: flow_fields
  type(bc_profile), allocatable :: profile

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

  ! Write gradients to solution file
  write_gradients = .true.

  ! Read and set BC profiles
  ! Read u componemt (1st column)
  call read_bc_profile(case_path // '.blasius.prf', 1, profile)
  profile%coordinates(:) = profile%coordinates(:) / mesh%geo%scalefactor
  profile%centre(:) = [ -4.0_ccs_real, 0.0_ccs_real, 0.5_ccs_real ] 
  
  ! Set to 3rd boundary condition (inlet)
  call set_bc_profile(u, profile, 3)

  ! Initialise velocity field
  if (irank == par_env%root) print *, "Initialise velocity field"
  call initialise_flow(par_env, run_options, flow_fields, get_init_flow, get_init_mass_flux)

  ! Solve using SIMPLE algorithm
  if (irank == par_env%root) print *, "Start SIMPLE"

  call activate_timestepping()
  call set_timestep(run_options%solve%dt)

  call run_solver(par_env, run_options, eval_sources, postproc_bfs, flow_fields)
  
  ! Clean-up

  call timer(end_time)

  call timer_print_all(par_env)
  call timer_export_csv(par_env)
  
  call nullify_mesh_object()
  ! Finalise MPI
  call cleanup_parallel_environment(par_env)

contains

  pure subroutine get_init_flow(loc_p, field_name, init_val)

    use types, only: cell_locator

    type(cell_locator), intent(in) :: loc_p
    character(len=*), intent(in) :: field_name
    real(ccs_real), intent(inout) :: init_val

    ! Silence compiler warnings
    associate(foo=>loc_p, bar=>field_name, baz=>init_val)
    end associate
    
  end subroutine get_init_flow
  
  pure subroutine get_init_mass_flux(loc_f, init_val)
    use types, only: face_locator
    type(face_locator), intent(in) :: loc_f
    real(ccs_real), intent(inout) :: init_val

    associate (foo => loc_f, bar => init_val)
    end associate

  end subroutine

  !> Subroutine to define case-specific postprocessing.
  subroutine postproc_bfs(par_env, flow_fields)

    class(parallel_environment), allocatable, intent(in) :: par_env
    type(fluid), intent(in) :: flow_fields

    ! All cases must define this, but if they don't require case-specific processing then simply
    ! make a no-op (use associate to silence unused variable compiler warnings)
    associate(foo => par_env, bar => flow_fields)
    end associate
    
  end subroutine postproc_bfs

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

end program bfs
