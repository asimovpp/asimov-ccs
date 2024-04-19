!v Module file core.mod
!
! The core module implements the core functionality to define and run a CCS case.

module core

  use kinds, only: ccs_int, ccs_real
  use types, only: fluid
  use parallel_types, only: parallel_environment

  implicit none

  private
  public :: initialise_mesh
  public :: run_solver
  public :: ccs_options
  integer(ccs_int), parameter, public :: read_input_mesh = 1
  integer(ccs_int), parameter, public :: build_mesh_2d = 2
  integer(ccs_int), parameter, public :: build_mesh_3d = 3

  !v Type to contain the configuration of the CCS run
  type :: ccs_options
    character(len=:), allocatable :: case_path
    character(len=:), allocatable :: case_name
    character(len=:), allocatable :: mesh_path
    integer(ccs_int) :: num_steps
    real(ccs_real) :: dt
    integer(ccs_int) :: write_frequency
    integer(ccs_int) :: it_start
    integer(ccs_int) :: it_end
    real(ccs_real) :: res_target
    integer(ccs_int) :: init_mesh_type
    integer(ccs_int) :: cps
    real(ccs_real) :: domain_size
  end type ccs_options


  interface

    !v Subroutine to initialise the mesh.
    !
    ! This is responsible for the building or reading of the mesh (whichever is specified in the config file)
    module subroutine initialise_mesh(par_env, shared_env, run_options)
      class(parallel_environment), intent(in), allocatable :: par_env
      class(parallel_environment), intent(in), allocatable :: shared_env
      type(ccs_options), intent(in) :: run_options
    end subroutine initialise_mesh

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
