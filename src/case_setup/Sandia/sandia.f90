!> Program file for Sandia case
program sandia
#include "ccs_macros.inc"

  use petscvec
  use petscsys

  use core
  use constants, only: ndim
  use meshing, only: nullify_mesh_object, get_local_num_cells
  use kinds, only: ccs_real, ccs_int, ccs_long
  use types, only: fluid
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer, &
                      is_root
  use parallel_types, only: parallel_environment
  use scalars, only: update_scalars
  use read_config, only: get_enable_cell_corrections, get_store_residuals
  use boundary_conditions, only: set_bc_profile
  use utils, only: str
  use fields, only: dealloc_fluid_fields
  use timers, only: timer_init, timer_register_start, timer_register, timer_start, timer_stop, &
                    timer_print, timer_get_time, timer_print_all, timer_export_csv

  implicit none

  class(parallel_environment), allocatable:: par_env
  class(parallel_environment), allocatable:: shared_env

  type(ccs_options) :: run_options

  integer(ccs_int):: irank  ! MPI rank ID
  integer(ccs_int):: isize  ! Size of MPI world

  integer(ccs_int):: timer_index_total
  integer(ccs_int):: timer_index_init

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

  if (is_root(par_env)) print *, "Starting ", run_options%paths%case_name, " case!"

  ! Read mesh from .geo file
  call initialise_mesh(par_env, shared_env, run_options)

  ! Initialise fields
  if (is_root(par_env)) print *, "Initialise fields"

  ! Initialise the fields
  call initialise_fields(par_env, run_options, flow_fields)

  ! XXX: coupling BCs could be built here

  ! Initialise velocity field
  if (is_root(par_env)) print *, "Initialise velocity field"
  call initialise_flow(par_env, run_options, flow_fields, get_init_flow, get_init_mass_flux)

  ! Solve using SIMPLE algorithm
  if (is_root(par_env)) print *, "Start SIMPLE"

  call timer_stop(timer_index_init)

  call run_solver(par_env, run_options, eval_sources, postproc_sandia, flow_fields)
  
  ! Clean-up

  call timer_stop(timer_index_total)

  call timer_print_all(par_env)
  call timer_export_csv(par_env)

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

    if (field_name == "scalar") then
      call get_centre(loc_p, x_p)
      if (x_p(1) < -0.08) then
        init_val = 1.0_ccs_real 
      else
        init_val = 0.0_ccs_real 
      end if
    else ! anything but scalar field
      init_val = init_val ! Accept whatever initial value is set
    end if

  end subroutine
  
  pure subroutine get_init_mass_flux(loc_f, init_val)
    use types, only: face_locator
    type(face_locator), intent(in) :: loc_f
    real(ccs_real), intent(inout) :: init_val

    associate (foo => loc_f, bar => init_val)

    end associate

  end subroutine

  !> Subroutine to define case-specific postprocessing.
  subroutine postproc_sandia(par_env, flow_fields)

    class(parallel_environment), allocatable, intent(in) :: par_env
    type(fluid), intent(in) :: flow_fields

    ! All cases must define this, but if they don't require case-specific processing then simply
    ! make a no-op (use associate to silence unused variable compiler warnings)
    associate(foo => par_env, bar => flow_fields)
    end associate
    
  end subroutine postproc_sandia

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

end program sandia
