!v Program file for LidDrivenCavity case
!
!  @build mpi+petsc

program ldc
#include "ccs_macros.inc"

  use petscvec
  use petscsys

  use core
  use case_config, only: write_gradients
  use constants, only: cell, face
  use kinds, only: ccs_real, ccs_int
  use types, only: fluid
  use parallel, only: initialise_parallel_environment, &
                      create_new_par_env, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments, sync, is_root
  use meshing, only: nullify_mesh_object
  use parallel_types, only: parallel_environment
  use vec, only: create_vector, set_vector_location, get_vector_data, restore_vector_data
  use timers, only: timer_init, timer_register_start, timer_register, timer_start, timer_stop, timer_print, timer_print_all
  use utils, only: dealloc_fluid_fields

  implicit none

  class(parallel_environment), allocatable :: par_env
  class(parallel_environment), allocatable :: shared_env

  integer(ccs_int) :: irank ! MPI rank ID
  integer(ccs_int) :: isize ! Size of MPI world

  integer(ccs_int) :: timer_index_init, timer_index_total

  type(fluid) :: flow_fields

  type(ccs_options) :: run_options

  ! Launch MPI
  call initialise_parallel_environment(par_env)
  call timer_init()

  irank = par_env%proc_id
  isize = par_env%num_procs

  call timer_register_start("Elapsed time", timer_index_total, is_total_time=.true.)
  call timer_register_start("Init time", timer_index_init)

  call get_config(par_env, run_options)
  call configure_parallelism(run_options, par_env, shared_env)
  if (irank == par_env%root) print *, "Starting ", run_options%paths%case_name, " case!"

  ! Create a mesh
  call initialise_mesh(par_env, shared_env, run_options)

  ! Initialise fields
  if (irank == par_env%root) print *, "Initialise fields"
  write_gradients = .false.
  call initialise_fields(par_env, run_options, flow_fields)

  ! Initialise velocity field
  if (irank == par_env%root) print *, "Initialise velocity field"
  call initialise_flow(par_env, run_options, flow_fields, get_init_flow, get_init_mass_flux)

  call timer_stop(timer_index_init)
  call run_solver(par_env, run_options, eval_sources, postproc_ldc, flow_fields)

  ! Clean-up
  call dealloc_fluid_fields(flow_fields)

  call timer_stop(timer_index_total)

  call timer_print_all(par_env)

  call nullify_mesh_object()

  ! Finalise MPI
  call cleanup_parallel_environment(par_env)

contains

  subroutine postproc_ldc(par_env, flow_fields)

    class(parallel_environment), allocatable, intent(in) :: par_env
    type(fluid), intent(in) :: flow_fields

    ! Silence compiler warnings
    associate(foo => par_env, bar => flow_fields)
    end associate
    
  end subroutine postproc_ldc

  pure subroutine get_init_flow(loc_p, field_name, init_val)

    use types, only: cell_locator
    use meshing, only: get_centre

    type(cell_locator), intent(in) :: loc_p
    character(len=*), intent(in) :: field_name
    real(ccs_real), intent(inout) :: init_val

    ! Silence ompiler warnings
    associate(foo => loc_p, bar => field_name, baz => init_val)
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
