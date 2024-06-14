!v Submodule implementing the core solver functionality.
!
! This provides the interface that runs the solver.

submodule (core) core_solver
#include "ccs_macros.inc"
  use kinds, only: ccs_int
  use types, only: fluid
  use parallel_types, only: parallel_environment
  use utils, only: debug_print

  use ccs_base, only: mesh
  
  use pv_coupling, only: solve_nonlinear

  use timers, only: timer_register, timer_start, timer_stop

  use timestepping, only: timestepping_is_active
  
  implicit none

  integer(ccs_int):: timer_index_sol
  integer(ccs_int):: timer_index_io_sol

contains

  module subroutine run_solver(par_env, run_options, postproc, flow_fields)

    class(parallel_environment), allocatable, intent(in) :: par_env
    type(ccs_options), intent(in) :: run_options
    interface
      subroutine postproc(par_env, flow_fields)
        use types, only: fluid
        use parallel_types, only: parallel_environment
        class(parallel_environment), allocatable, intent(in) :: par_env
        type(fluid), intent(in) :: flow_fields
      end subroutine
    end interface
    type(fluid), intent(inout) :: flow_fields

    integer(ccs_int) :: t ! Timestep counter
    integer(ccs_int) :: num_steps

    integer(ccs_int) :: it_start, it_end

    logical:: diverged = .false.

    integer(ccs_int) :: write_frequency
    
    it_start = run_options%solve%it_start
    it_end = run_options%solve%it_end
    if (timestepping_is_active()) then
      num_steps = run_options%solve%num_steps
    else
      ! num steps may not have been set
      num_steps = 1
    end if
    
    call timer_register("I/O time for solution", timer_index_io_sol)
    call timer_register("Solver time inc I/O", timer_index_sol)
    
    write_frequency = run_options%io%write_frequency
    do t = 1, num_steps
      call timer_start(timer_index_sol)

      ! XXX: Coupler update here
      
      call solve_nonlinear(par_env, run_options, mesh, flow_fields, diverged)

      ! XXX: Or coupler update here?
      
      if (par_env%proc_id == par_env%root) then
        print *, "TIME = ", t
      end if

      call postproc(par_env, flow_fields)

      ! If a STOP file exist, write solution and exit the main simulation loop
      if (check_stop_run(par_env, diverged)) then
        call write_step(par_env, run_options, t, flow_fields)
        call dprint("STOP run")
        exit
      end if

      if (check_to_write(run_options, t)) then
        call write_step(par_env, run_options, t, flow_fields)
      end if
      call timer_stop(timer_index_sol)
    end do

  end subroutine run_solver

  logical function check_stop_run(par_env, diverged)

    use parallel, only: query_stop_run
    
    class(parallel_environment), intent(in) :: par_env
    logical, intent(in) :: diverged
    
    check_stop_run = (query_stop_run(par_env) .or. diverged)

  end function check_stop_run
  
  logical pure function check_to_write(run_options, t)

    type(ccs_options), intent(in) :: run_options
    integer(ccs_int), intent(in) :: t
    
    associate(num_steps => run_options%solve%num_steps, &
              write_frequency => run_options%io%write_frequency)
      check_to_write = ((t == 1) .or. (t == num_steps) .or. (mod(t, write_frequency) == 0))
    end associate
    
  end function check_to_write
  
  subroutine write_step(par_env, run_options, t, flow_fields)

    use io_visualisation, only: write_solution

    class(parallel_environment), allocatable, intent(in) :: par_env
    type(ccs_options), intent(in) :: run_options
    integer(ccs_int), intent(in) :: t
    type(fluid), intent(inout) :: flow_fields

    character(len=:), allocatable :: case_path
    integer(ccs_int) :: num_steps
    real(ccs_real) :: dt
    
    case_path = run_options%paths%case_path
    num_steps = run_options%solve%num_steps
    dt = run_options%solve%dt
    
    call timer_start(timer_index_io_sol)
    call write_solution(par_env, case_path, mesh, flow_fields, t, num_steps, dt)
    call timer_stop(timer_index_io_sol)

  end subroutine write_step

end submodule core_solver

