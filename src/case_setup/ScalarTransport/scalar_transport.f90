!v Program file for ScalarTransport case
!
!  @build mpi+petsc

program scalar_transport
#include "ccs_macros.inc"

  use petscvec
  use petscsys

  use core
  use kinds, only: ccs_real, ccs_int
  use types, only: field, ccs_mesh, fluid
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer, &
                      is_root
  use parallel_types, only: parallel_environment
  use mesh_utils, only: build_mesh
  use meshing, only: get_global_num_cells, get_centre, count_neighbours, &
                     create_cell_locator, create_face_locator, create_neighbour_locator, &
                     get_local_index, get_boundary_status, get_face_normal, nullify_mesh_object, &
                     get_local_num_cells
  use scalars, only: update_scalars
  use utils, only: exit_print, add_field_to_outputlist
  use fields, only: get_field, dealloc_fluid_fields
  use profiler, only: profiler_init, profiler_shutdown, profiler_begin_region, profiler_end_region

  implicit none

  class(parallel_environment), allocatable :: par_env
  class(parallel_environment), allocatable, target :: shared_env

  type(fluid) :: flow_fields

  type(ccs_options) :: run_options
  
  ! Launch MPI
  call initialise_parallel_environment(par_env)
  call profiler_init()

  call get_config(par_env, run_options)
  call configure_parallelism(run_options, par_env, shared_env)

  call profiler_begin_region('Total elapsed time')

  call profiler_begin_region('Total initialisation')

  if (is_root(par_env)) print *, "Starting ", run_options%paths%case_name, " case!"

  call initialise_mesh(par_env, shared_env, run_options)

  ! Initialise fields
  if (is_root(par_env)) print *, "Initialise fields"

  call initialise_fields(par_env, run_options, flow_fields)

  ! Initialise velocity field
  if (is_root(par_env)) print *, "Initialise flow field"
  call initialise_flow(par_env, run_options, flow_fields, get_init_flow, get_init_mass_flux)

  ! Solve using SIMPLE algorithm
  if (is_root(par_env)) print *, "Start scalar solver"

  call profiler_end_region('Total initialisation')

  ! Write out mesh and solution
  call run_solver(par_env, run_options, eval_sources, postproc_scalar, flow_fields)

  call profiler_end_region('Total elapsed time')
  call profiler_shutdown(par_env)

  ! Clean-up
  call dealloc_fluid_fields(flow_fields)

  ! Finalise MPI
  call nullify_mesh_object()
  call cleanup_parallel_environment(par_env)

contains

  pure subroutine get_init_flow(loc_p, field_name, init_val)

    use types, only: cell_locator

    type(cell_locator), intent(in) :: loc_p
    character(len=*), intent(in) :: field_name
    real(ccs_real), intent(inout) :: init_val

    real(ccs_real), dimension(3) :: x
    real(ccs_real), dimension(2) :: r

    real(ccs_real) :: whisky

    real(ccs_real) :: L, c

    L = run_options%mesh%domain_size
    c = L / 2

    call get_centre(loc_p, x)
    r = x(1:2) - c

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

    real(ccs_real), dimension(3) :: face_normal, v
    real(ccs_real) :: rmag
    real(ccs_real) :: theta
    logical :: is_boundary

    real(ccs_real), dimension(3) :: x
    real(ccs_real), dimension(2) :: r
    real(ccs_real) :: L, c
    
    L = run_options%mesh%domain_size
    c = L / 2

    call get_boundary_status(loc_f, is_boundary)
    if (is_boundary) then
      init_val = 0.0_ccs_real
    else
      call get_face_normal(loc_f, face_normal)
      call get_centre(loc_f, x)

      ! Create a rotating field that decays to zero on the boundaries
      r = x(1:2) - c
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

end program scalar_transport
