!v Module file core.mod
!
! The core module implements the core functionality to define and run a CCS case.

module core

  use kinds, only: ccs_int, ccs_real
  use types, only: fluid
  use parallel_types, only: parallel_environment

  implicit none

  private

  public :: run_solver
  public :: ccs_options

  !v Type to contain the configuration of the CCS run
  type :: ccs_options
    character(len=:), allocatable :: case_path
    integer(ccs_int) :: num_steps
    real(ccs_real) :: dt
    integer(ccs_int) :: write_frequency
    integer(ccs_int) :: it_start
    integer(ccs_int) :: it_end
    real(ccs_real) :: res_target
  end type ccs_options

  interface
    !v Subroutine to run a flow problem.
    !
    ! This is responsible for the time loop, calling post-processing subroutines and performing
    ! solution output.
    module subroutine run_solver(par_env, run_options, flow_fields)
      class(parallel_environment), allocatable, intent(in) :: par_env
      type(ccs_options), intent(in) :: run_options
      type(fluid), intent(inout) :: flow_fields
    end subroutine run_solver
  end interface
  
end module core
