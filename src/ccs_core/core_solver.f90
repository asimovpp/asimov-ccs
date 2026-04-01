!v Submodule implementing the core solver functionality.
!
! This provides the interface that runs the solver.

submodule (core) core_solver
#include "ccs_macros.inc"
  use kinds, only: ccs_int
  use types, only: fluid
  use parallel, only: is_root
  use parallel_types, only: parallel_environment
  use utils, only: debug_print

  use ccs_base, only: mesh
  
  use pv_coupling, only: solve_nonlinear
  use scalars, only: update_scalars
  use profiler, only: profiler_begin_region, profiler_end_region
  use timestepping, only: timestepping_is_active, finalise_timestep
  
  implicit none

contains

  !v Subroutine to run a flow problem.
  !
  ! This is responsible for the time loop, calling post-processing subroutines and performing
  ! solution output.
  module subroutine run_solver(par_env, run_options, eval_sources, postproc, flow_fields)

    use timestepping, only: activate_timestepping, set_timestep
    use flow_stats, only: report_cfl

    class(parallel_environment), allocatable, intent(in) :: par_env !< The parallel environment
    type(ccs_options), intent(in) :: run_options                    !< The runtime configuration
    interface
      !v Subroutine to evaluate source terms, case-specific.
      !
      !  Note this should return the integrated source.
      subroutine eval_sources(flow, phi, R, S)
        use types, only: fluid, field, ccs_vector
        type(fluid), intent(in) :: flow !< Provides access to full flow field
        class(field), intent(in) :: phi !< Field being transported
        class(ccs_vector), intent(inout) :: R !< Work vector (for evaluating linear/implicit sources)
        class(ccs_vector), intent(inout) :: S !< Work vector (for evaluating fixed/explicit sources)
      end subroutine eval_sources
    end interface
    interface
      !v Subroutine to perform online analysis of the solution, case-specific.
      subroutine postproc(par_env, flow_fields)
        use types, only: fluid
        use parallel_types, only: parallel_environment
        class(parallel_environment), allocatable, intent(in) :: par_env !< The parallel environment
        type(fluid), intent(in) :: flow_fields                          !< The flow field structure
      end subroutine
    end interface
    type(fluid), intent(inout) :: flow_fields !< The flow field structure, contains the solution

    integer(ccs_int) :: t ! Timestep counter
    integer(ccs_int) :: num_steps

    integer(ccs_int) :: it_start, it_end

    logical:: diverged = .false.
    
    if (run_options%solve%unsteady) then
      call activate_timestepping()
      call set_timestep(run_options%solve%dt)
    end if
  
    it_start = run_options%solve%it_start
    it_end = run_options%solve%it_end
    if (timestepping_is_active()) then
      num_steps = run_options%solve%num_steps
    else
      ! num steps may not have been set
      num_steps = 1
    end if
    
    do t = 1, num_steps
      call profiler_begin_region("Solver time inc I/O")

      ! XXX: Coupler update here
      call report_cfl(par_env, flow_fields)
      call advance_step(par_env, run_options, eval_sources, flow_fields, diverged)
      ! XXX: Or coupler update here?
      
      if (timestepping_is_active()) then
        if (is_root(par_env)) then
          print *, "TIME = ", t
        end if
      end if
      
      call postproc(par_env, flow_fields)

      if (check_stop_run(par_env, run_options, t, flow_fields, diverged)) then
        call profiler_end_region("Solver time inc I/O")
        exit
      end if

      if (check_to_write(run_options, t)) then
        call write_step(par_env, run_options, t, flow_fields)
      end if

      call profiler_end_region("Solver time inc I/O")
    end do


  end subroutine run_solver

  !> Advance the simulation by one step.
  subroutine advance_step(par_env, run_options, eval_sources, flow_fields, diverged)

    class(parallel_environment), allocatable, intent(in) :: par_env
    type(ccs_options), intent(in) :: run_options
    interface
      !v Subroutine to evaluate source terms, case-specific.
      !
      !  Note this should return the integrated source.
      subroutine eval_sources(flow, phi, R, S)
        use types, only: fluid, field, ccs_vector
        type(fluid), intent(in) :: flow !< Provides access to full flow field
        class(field), intent(in) :: phi !< Field being transported
        class(ccs_vector), intent(inout) :: R !< Work vector (for evaluating linear/implicit sources)
        class(ccs_vector), intent(inout) :: S !< Work vector (for evaluating fixed/explicit sources)
      end subroutine eval_sources
    end interface
    type(fluid), intent(inout) :: flow_fields
    logical, intent(out) :: diverged

    logical :: flow_sol

    flow_sol = check_flow_sol(par_env, flow_fields)
      
    if (flow_sol) then
      call solve_nonlinear(par_env, run_options, eval_sources, mesh, flow_fields, diverged)
    else
      ! Only scalar transport
      diverged = .false.
      call update_scalars(par_env, mesh, eval_sources, flow_fields)
    end if
    call finalise_timestep()
    
  end subroutine advance_step

  !> Checks for stop conditions.
  logical function check_stop_run(par_env, run_options, t, flow_fields, diverged)

    class(parallel_environment), intent(in), allocatable :: par_env
    type(ccs_options), intent(in) :: run_options
    integer, intent(in) :: t
    type(fluid), intent(inout) :: flow_fields
    logical, intent(in) :: diverged

    if (stop_if_diverged(par_env, run_options, t, flow_fields, diverged)) then
      check_stop_run = .true.
    else if (stop_on_request(par_env, run_options, t, flow_fields)) then
      check_stop_run = .true.
    else
      check_stop_run = .false.
    end if
    
  end function check_stop_run
  
  !> Checks for stop condition if divergence is detected and dumps solution if this occurs.
  logical function stop_if_diverged(par_env, run_options, t, flow_fields, diverged)

    class(parallel_environment), intent(in), allocatable :: par_env
    type(ccs_options), intent(in) :: run_options
    integer, intent(in) :: t
    type(fluid), intent(inout) :: flow_fields
    logical, intent(in) :: diverged
    
    if (diverged) then
      if (is_root(par_env)) then
        print *, "INFO: Divergence detected"
      end if
      call dump_run(par_env, run_options, t, flow_fields)

      stop_if_diverged = .true.
    else
      stop_if_diverged = .false.
    end if
    
  end function stop_if_diverged

  !> Checks for stop condition due to STOP file and dumps solution if this occurs.
  logical function stop_on_request(par_env, run_options, t, flow_fields)

    use parallel, only: query_stop_run

    class(parallel_environment), intent(in), allocatable :: par_env
    type(ccs_options), intent(in) :: run_options
    integer, intent(in) :: t
    type(fluid), intent(inout) :: flow_fields

    if (query_stop_run(par_env)) then
      if (is_root(par_env)) then
        print *, "INFO: Found STOP file"
      end if
      call dump_run(par_env, run_options, t, flow_fields)

      stop_on_request = .true.
    else
      stop_on_request = .false.
    end if

  end function stop_on_request

  !> Dumps the solution when stopping a simulation.
  subroutine dump_run(par_env, run_options, t, flow_fields)

    class(parallel_environment), intent(in), allocatable :: par_env
    type(ccs_options), intent(in) :: run_options
    integer, intent(in) :: t
    type(fluid), intent(inout) :: flow_fields

    if (is_root(par_env)) then
      print *, "STOPPING SIMULATION"
    end if
    call write_step(par_env, run_options, t, flow_fields)

  end subroutine dump_run

  !v Check whether we are solving fluid flow or scalars only. If pressure and at least one of u,v,w
  !  are present then we are solving the flow field, otherwise it is scalar transport with frozen
  !  flow field.
  logical function check_flow_sol(par_env, flow_fields) 

    use types, only: field
    use parallel, only: is_root
    use fields, only: get_field
    
    class(parallel_environment), intent(in) :: par_env !< The parallel environemnt
    type(fluid), intent(in) :: flow_fields             !< The flow field structure

    logical :: have_p, have_vel
    integer :: i_field
    class(field), pointer :: phi

    have_p = .false.
    have_vel = .false.

    do i_field = 1, size(flow_fields%fields)
      call get_field(flow_fields, i_field, phi)
      if (phi%name == "p" .and. phi%solver_parameters%solve) then
        have_p = .true.
      else if (((phi%name == "u") .or. (phi%name == "v") .or. (phi%name == "w")) .and. phi%solver_parameters%solve) then
        have_vel = .true.
      end if
    end do

    check_flow_sol = have_p .and. have_vel
    if (is_root(par_env)) then
      if (check_flow_sol) then
        print *, "Solving fluid flow"
      else
        print *, "Solving scalar transport only"
      end if
    end if
    
  end function check_flow_sol
  
  !> Predicate to test if conditions for solution output are met
  logical pure function check_to_write(run_options, t)

    type(ccs_options), intent(in) :: run_options !< The runtime configuration
    integer(ccs_int), intent(in) :: t            !< The timestep counter
    
    if (timestepping_is_active()) then
      associate(num_steps => run_options%solve%num_steps, &
                write_frequency => run_options%io%write_frequency)
        check_to_write = ((t == 1) .or. (t == num_steps) .or. (mod(t, write_frequency) == 0))
      end associate
    else
      ! End of steady run
      check_to_write = .true.
    end if
    
  end function check_to_write
  
  !> Utility subroutine to write the solution for a step
  subroutine write_step(par_env, run_options, t, flow_fields)

    use io_visualisation, only: write_solution

    class(parallel_environment), allocatable, intent(in) :: par_env !< The parallel environment
    type(ccs_options), intent(in) :: run_options                    !< The runtime configuration
    integer(ccs_int), intent(in) :: t                               !< The timestep counter
    type(fluid), intent(inout) :: flow_fields                       !< The flow field structure

    integer(ccs_int) :: num_steps
    real(ccs_real) :: dt
    
    num_steps = run_options%solve%num_steps
    dt = run_options%solve%dt
    
    call profiler_begin_region("I/O time for solution")
    if (timestepping_is_active()) then
      call write_solution(par_env, run_options, mesh, flow_fields, t, num_steps, dt)
    else
      call write_solution(par_env, run_options, mesh, flow_fields)
    end if
    call profiler_end_region("I/O time for solution")

  end subroutine write_step

end submodule core_solver

