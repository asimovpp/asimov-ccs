!> Controlled solver driver for boundary mass-flux reporting tests.
!> The LIT driver runs four Input/BoundaryMassFlux*_config.yaml configurations on four ranks:
!> Enabled checks one root-only report after each of two steps, at times 0.1 and 0.2;
!> Steady checks one report after the solve; Disabled and Absent check that no report is printed.
!> Boundary flux starts at one and is doubled by post-processing after each report. On the unit
!> square this gives patch integrals of one then two, and global imbalances of four then eight,
!> verifying that each report uses the current flux. The flag controls output, not integration.
program boundary_mass_flux_reporting_driver
#include "ccs_macros.inc"

  use core
  use fields, only: dealloc_fluid_fields
  use kinds, only: ccs_real
  use meshing, only: get_boundary_status, nullify_mesh_object
  use parallel, only: cleanup_parallel_environment, initialise_parallel_environment
  use parallel_types, only: parallel_environment
  use profiler, only: profiler_begin_region, profiler_end_region, profiler_init, profiler_shutdown
  use types, only: field, fluid

  implicit none

  class(parallel_environment), allocatable :: par_env
  class(parallel_environment), allocatable, target :: shared_env
  type(fluid) :: flow_fields
  type(ccs_options) :: run_options

  call initialise_parallel_environment(par_env)
  call profiler_init()
  call get_config(par_env, run_options)
  call configure_parallelism(run_options, par_env, shared_env)

  call profiler_begin_region("Total elapsed time")
  call profiler_begin_region("Total initialisation")
  call initialise_mesh(par_env, shared_env, run_options)
  call initialise_fields(par_env, run_options, flow_fields)
  call initialise_flow(par_env, run_options, flow_fields, get_init_flow, &
                       get_init_mass_flux)
  call profiler_end_region("Total initialisation")

  call run_solver(par_env, run_options, eval_sources, scale_mass_flux, flow_fields)

  call profiler_end_region("Total elapsed time")
  call profiler_shutdown(par_env)
  call dealloc_fluid_fields(flow_fields)
  call nullify_mesh_object()
  call cleanup_parallel_environment(par_env)

contains

  pure subroutine get_init_flow(loc_p, field_name, init_val)
    use types, only: cell_locator

    type(cell_locator), intent(in) :: loc_p
    character(len=*), intent(in) :: field_name
    real(ccs_real), intent(inout) :: init_val

    associate (unused_loc => loc_p, unused_name => field_name)
      init_val = 0.0_ccs_real
    end associate
  end subroutine get_init_flow

  pure subroutine get_init_mass_flux(loc_f, init_val)
    use types, only: face_locator

    type(face_locator), intent(in) :: loc_f
    real(ccs_real), intent(inout) :: init_val

    logical :: is_boundary

    call get_boundary_status(loc_f, is_boundary)
    init_val = merge(1.0_ccs_real, 0.0_ccs_real, is_boundary)
  end subroutine get_init_mass_flux

  subroutine scale_mass_flux(par_env, flow_fields)
    use fields, only: get_field
    use utils, only: update
    use vec, only: get_vector_data, restore_vector_data

    class(parallel_environment), allocatable, intent(in) :: par_env
    type(fluid), intent(in) :: flow_fields

    class(field), pointer :: mass_flux
    real(ccs_real), dimension(:), pointer :: mass_flux_values

    call get_field(flow_fields, "mf", mass_flux)
    call get_vector_data(mass_flux%values, mass_flux_values)
    mass_flux_values = 2.0_ccs_real * mass_flux_values
    call restore_vector_data(mass_flux%values, mass_flux_values)
    call update(mass_flux%values)

    associate (unused_env => par_env)
    end associate
  end subroutine scale_mass_flux

  subroutine eval_sources(flow, phi, R, S)
    use fv, only: zero_sources
    use types, only: ccs_vector

    type(fluid), intent(in) :: flow
    class(field), intent(in) :: phi
    class(ccs_vector), intent(inout) :: R
    class(ccs_vector), intent(inout) :: S

    call zero_sources(flow, phi, R, S)
  end subroutine eval_sources

end program boundary_mass_flux_reporting_driver
