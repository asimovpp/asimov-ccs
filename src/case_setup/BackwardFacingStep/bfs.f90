!> Program file for BackwardFacingStep case
program bfs
#include "ccs_macros.inc"

  use petscvec
  use petscsys

  use core
  use ccs_base, only: mesh
  use kinds, only: ccs_real, ccs_int, ccs_long
  use types, only: field, fluid, bc_profile
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer, &
                      is_root
  use parallel_types, only: parallel_environment
  use fields, only: get_field, dealloc_fluid_fields
  use boundary_conditions, only: set_bc_profile
  use meshing, only: nullify_mesh_object
  use utils, only: str
  use profiler, only: profiler_init, profiler_shutdown, profiler_begin_region, profiler_end_region
  use logging, only: log_unit_out

  implicit none

  class(parallel_environment), allocatable :: par_env
  class(parallel_environment), allocatable :: shared_env

  class(field), pointer :: u

  type(fluid) :: flow_fields
  type(bc_profile), allocatable :: profile

  type(ccs_options) :: run_options

  ! Launch MPI
  call initialise_parallel_environment(par_env)
  call profiler_init()

  call get_config(par_env, run_options)
  call configure_parallelism(run_options, par_env, shared_env)

  call profiler_begin_region('Total elapsed time')

  call profiler_begin_region('Total initialisation')

  if (is_root(par_env)) write(log_unit_out,*) "Starting ", run_options%paths%case_name, " case!"

  call initialise_mesh(par_env, shared_env, run_options)

  ! Initialise fields
  if (is_root(par_env)) write(log_unit_out,*) "Initialise fields"
  call initialise_fields(par_env, run_options, flow_fields)

  ! Read and set BC profiles
  ! Read u componemt (1st column)
  call get_field(flow_fields, "u", u)
  call read_bc_profile(run_options%paths%case_path // '.blasius.prf', 1, profile)
  profile%coordinates(:) = profile%coordinates(:) / mesh%geo%scalefactor
  profile%centre(:) = [ -4.0_ccs_real, 0.0_ccs_real, 0.5_ccs_real ] 
  
  ! Set to 3rd boundary condition (inlet)
  call set_bc_profile(u, profile, 3)
  nullify(u)

  ! Initialise velocity field
  if (is_root(par_env)) write(log_unit_out,*) "Initialise velocity field"
  call initialise_flow(par_env, run_options, flow_fields, get_init_flow, get_init_mass_flux)

  ! Solve using SIMPLE algorithm
  call profiler_end_region('Total initialisation')
  call run_solver(par_env, run_options, eval_sources, postproc_bfs, flow_fields)
  call profiler_end_region('Total elapsed time')

  call profiler_shutdown(par_env)
  
  call dealloc_fluid_fields(flow_fields)
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
    character(len = 128) :: header_string, tmp
    integer(ccs_int) :: num_field, i
    integer :: io_err, unit_io

    allocate(profile)

    allocate(profile%centre(3))
    allocate(profile%values(0))
    allocate(profile%coordinates(0))

    open(newunit = unit_io, file = trim(filename), status='old', action='read')

    read(unit_io, *)                      ! ignore profile type
    read(unit_io, *) tmp, profile%centre  ! read centre
    read(unit_io, *)                      ! ignore tolerance
    read(unit_io, *)                      ! ignore scaling
    read(unit_io, '(A)') header_string

    ! Count the number of fields in file
    num_field = -1
    do i = 1, len(header_string)
      if (header_string(i:i) == ',') then
        num_field = num_field + 1
      end if
    end do

    allocate(tmp_values(num_field))

    ! Read file profile table
    do while (.true.)

      read(unit_io, *, iostat = io_err) tmp_coord, tmp_values
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
