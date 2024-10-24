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

  use timers, only: timer_register, timer_start, timer_stop

  use timestepping, only: timestepping_is_active
  
  implicit none

  integer(ccs_int):: timer_index_sol
  integer(ccs_int):: timer_index_io_sol

contains

  module subroutine run_solver(par_env, run_options, eval_sources, postproc, flow_fields)

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

    logical :: flow_sol
    
    it_start = run_options%solve%it_start
    it_end = run_options%solve%it_end
    if (timestepping_is_active()) then
      num_steps = run_options%solve%num_steps
    else
      ! num steps may not have been set
      num_steps = 1
    end if

    flow_sol = check_flow_sol(par_env, flow_fields)
    
    call timer_register("I/O time for solution", timer_index_io_sol)
    call timer_register("Solver time inc I/O", timer_index_sol)
    
    write_frequency = run_options%io%write_frequency
    do t = 1, num_steps
      call timer_start(timer_index_sol)

      ! XXX: Coupler update here
      
      if (flow_sol) then
        call solve_nonlinear(par_env, run_options, eval_sources, mesh, flow_fields, diverged)
      else
        ! Only scalar transport
        call update_scalars(par_env, mesh, eval_sources, flow_fields)
      end if

      ! XXX: Or coupler update here?
      
      if (timestepping_is_active()) then
        if (is_root(par_env)) then
          print *, "TIME = ", t
        end if
      end if
      
      call postproc(par_env, flow_fields)

      ! If a STOP file exist, write solution and exit the main simulation loop
      if (check_stop_run(par_env, diverged)) then
        if (is_root(par_env)) then
          if (diverged) then
            print *, "INFO: Divergence detected, stopping"
          else
            print *, "INFO: STOP file present, stopping"
          end if
        end if
        call write_step(par_env, run_options, t, flow_fields)
        exit
      end if

      if (check_to_write(run_options, t)) then
        call write_step(par_env, run_options, t, flow_fields)
      end if
      call timer_stop(timer_index_sol)
    end do

  end subroutine run_solver

  !v Check whether we are solving fluid flow or scalars only. If pressure and at least one of u,v,w
  !  are present then we are solving the flow field, otherwise it is scalar transport with frozen
  !  flow field.
  logical function check_flow_sol(par_env, flow_fields) 

    use types, only: field
    use parallel, only: is_root
    use utils, only: get_field
    
    class(parallel_environment), intent(in) :: par_env
    type(fluid), intent(in) :: flow_fields

    logical :: have_p, have_vel
    integer :: i_field
    class(field), pointer :: phi

    have_p = .false.
    have_vel = .false.

    do i_field = 1, size(flow_fields%fields)
      call get_field(flow_fields, i_field, phi)
      if (phi%name == "p") then
        have_p = .true.
      else if ((phi%name == "u") .or. (phi%name == "v") .or. (phi%name == "w")) then
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

  logical function check_stop_run(par_env, diverged)

    use parallel, only: query_stop_run
    
    class(parallel_environment), intent(in) :: par_env
    logical, intent(in) :: diverged
    
    check_stop_run = (query_stop_run(par_env) .or. diverged)

  end function check_stop_run
  
  logical pure function check_to_write(run_options, t)

    type(ccs_options), intent(in) :: run_options
    integer(ccs_int), intent(in) :: t
    
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
    if (timestepping_is_active()) then
      call write_solution(par_env, run_options, mesh, flow_fields, t, num_steps, dt)
    else
      call write_solution(par_env, run_options, mesh, flow_fields)
    end if
    call timer_stop(timer_index_io_sol)

  end subroutine write_step

end submodule core_solver

