!> @brief Regression test for restart target-step semantics
!
!> @description Verifies that, on a restart, the configured `steps` value is the final
!>              timestep to run up to (an absolute step number) rather than an additional
!>              number of steps on top of the restart step. The test first runs the
!>              Taylor-Green vortex case for `restart_step` timesteps (writing the
!>              solution), then restarts from that step with `steps = target_step` and
!>              checks that the solver stops at `target_step` (not `target_step + restart_step`).
program test_restart_target
#include "ccs_macros.inc"

  use core
  use testing_lib
  use tgv2d_core, only: get_init_flow, get_init_mass_flux, eval_sources, postproc_tgv
  use timestepping, only: get_current_step, reset_timestepping
  use fields, only: dealloc_fluid_fields
  use io_visualisation, only: reset_io_visualisation
  use utils, only: str, reset_outputlist_counter
  use types, only: fluid
  use kinds, only: ccs_int

  implicit none

  integer(ccs_int), parameter :: restart_step = 10_ccs_int  ! step to resume from
  integer(ccs_int), parameter :: target_step = 20_ccs_int   ! final step to run up to
  integer(ccs_int), parameter :: cps = 8_ccs_int            ! small mesh for a fast test

  type(fluid) :: flow_fields
  type(ccs_options) :: run_options
  integer(ccs_int) :: final_step
  character(len=:), allocatable :: restart_file

  call init()

  call get_config(par_env, run_options)

  ! Build the mesh once; it is reused for both phases.
  run_options%mesh%cps = cps
  call initialise_mesh(par_env, shared_env, run_options)

  ! ---- Phase 1: run restart_step timesteps (no restart), writing the solution ----
  run_options%variables%restart = .false.
  run_options%solve%num_steps = restart_step
  call initialise_fields(par_env, run_options, flow_fields)
  call initialise_flow(par_env, run_options, flow_fields, get_init_flow, get_init_mass_flux)
  call run_solver(par_env, run_options, eval_sources, postproc_tgv, flow_fields)

  ! Reset to a clean state before the restart phase.
  call reset_timestepping()
  call reset_outputlist_counter()
  call reset_io_visualisation()
  call dealloc_fluid_fields(flow_fields)

  ! ---- Phase 2: restart from restart_step and run up to target_step ----
  restart_file = run_options%paths%case_name // "_sol_" // str(restart_step) // ".h5"
  run_options%variables%restart = .true.
  run_options%variables%restart_file = restart_file
  run_options%variables%restart_step = restart_step
  run_options%solve%num_steps = target_step
  call initialise_fields(par_env, run_options, flow_fields)
  call initialise_flow(par_env, run_options, flow_fields, get_init_flow, get_init_mass_flux)
  call run_solver(par_env, run_options, eval_sources, postproc_tgv, flow_fields)

  call get_current_step(final_step)
  call assert_eq(final_step, target_step, &
                 "On restart, 'steps' must be the final step to run up to (no extra steps beyond the target)")

  ! Clean up.
  call reset_timestepping()
  call reset_outputlist_counter()
  call reset_io_visualisation()
  call dealloc_fluid_fields(flow_fields)
  call finalise_mesh(par_env, .true.)

  call fin()

end program test_restart_target
