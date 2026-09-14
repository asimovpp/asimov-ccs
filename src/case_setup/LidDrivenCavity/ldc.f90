!v Program file for LidDrivenCavity case
!
!  @build mpi+petsc

program ldc
#include "ccs_macros.inc"

  use core
  use kinds, only: ccs_real
  use types, only: fluid
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, &
                      is_root
  use parallel_types, only: parallel_environment
  use fields, only: dealloc_fluid_fields
  use profiler, only: profiler_init, profiler_shutdown, profiler_begin_region, profiler_end_region
  use logging, only: log_unit_out

  implicit none

  class(parallel_environment), allocatable :: par_env
  class(parallel_environment), allocatable :: shared_env

  type(fluid) :: flow_fields

  type(ccs_options) :: run_options

  ! Launch MPI
  call initialise_parallel_environment(par_env)
  call profiler_init()

  call profiler_begin_region('Total elapsed time')
  call profiler_begin_region('Total initialisation')

  call get_config(par_env, run_options)
  call configure_parallelism(run_options, par_env, shared_env)
  if (is_root(par_env)) write (log_unit_out, *) "Starting ", run_options%paths%case_name, " case!"

  ! Create a mesh
  call initialise_mesh(par_env, shared_env, run_options)

  ! Initialise fields
  if (is_root(par_env)) write (log_unit_out, *) "Initialise fields"
  call initialise_fields(par_env, run_options, flow_fields)

  ! Initialise velocity field
  if (is_root(par_env)) write (log_unit_out, *) "Initialise velocity field"
  call initialise_flow(par_env, run_options, flow_fields, get_init_flow, get_init_mass_flux)

  call profiler_end_region('Total initialisation')
  call run_solver(par_env, run_options, eval_sources, postproc_ldc, flow_fields)

  ! Clean-up
  call dealloc_fluid_fields(flow_fields)

  call profiler_end_region('Total elapsed time')
  call profiler_shutdown(par_env)

  call finalise_mesh(par_env)

  ! Finalise MPI
  call cleanup_parallel_environment(par_env)

contains

  subroutine postproc_ldc(par_env, flow_fields)

    class(parallel_environment), allocatable, intent(in) :: par_env
    type(fluid), intent(in) :: flow_fields

    ! Silence compiler warnings
    associate (foo => par_env, bar => flow_fields)
    end associate

  end subroutine postproc_ldc

  pure subroutine get_init_flow(loc_p, field_name, init_val)

    use types, only: cell_locator

    type(cell_locator), intent(in) :: loc_p
    character(len=*), intent(in) :: field_name
    real(ccs_real), intent(inout) :: init_val

    ! Silence ompiler warnings
    associate (foo => loc_p, bar => field_name, baz => init_val)
    end associate

  end subroutine get_init_flow

  pure subroutine get_init_mass_flux(loc_f, init_val)

    use types, only: face_locator
    type(face_locator), intent(in) :: loc_f
    real(ccs_real), intent(inout) :: init_val

    ! Silence compiler warnings
    associate (foo => loc_f, bar => init_val)
    end associate

  end subroutine get_init_mass_flux

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

end program ldc
