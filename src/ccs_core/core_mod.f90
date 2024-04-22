!v Module file core.mod
!
! The core module implements the core functionality to define and run a CCS case.

module core

  use constants, only: ccs_string_len
  use kinds, only: ccs_int, ccs_real
  use types, only: fluid
  use parallel_types, only: parallel_environment

  implicit none

  private

  public :: get_config
  public :: configure_parallelism
  public :: initialise_mesh
  public :: run_solver
  public :: ccs_options

  integer(ccs_int), parameter, public :: read_input_mesh = 1
  integer(ccs_int), parameter, public :: build_mesh_2d = 2
  integer(ccs_int), parameter, public :: build_mesh_3d = 3

  !v Type to contain path-related options for the CCS run
  type :: ccs_paths
    character(len=:), allocatable :: case_name
    character(len=:), allocatable :: case_path
    character(len=:), allocatable :: mesh_path
    character(len=:), allocatable :: ccs_config_file
  end type ccs_paths
  
  !v Options for the mesh configuration
  type :: mesh_options
    integer(ccs_int) :: init_mesh_type
    integer(ccs_int) :: cps
    real(ccs_real) :: domain_size
  end type mesh_options
  
  !v Options for IO configuration
  type :: io_options
    integer(ccs_int) :: write_frequency
  end type io_options

  !v Options for variable declarations
  type :: variable_options
    character(len=ccs_string_len), dimension(:), allocatable :: variable_names
    integer(ccs_int), dimension(:), allocatable :: variable_types
  end type variable_options

  !v Options for solver configuration
  type :: solver_options
    integer(ccs_int) :: num_steps
    integer(ccs_int) :: num_iters
    integer(ccs_int) :: it_start
    integer(ccs_int) :: it_end
    real(ccs_real) :: dt
    real(ccs_real) :: res_target
    real(ccs_real) :: velocity_relax
    real(ccs_real) :: pressure_relax
  end type solver_options
  
  !v Type to contain the configuration of the CCS run
  type :: ccs_options
    type(ccs_paths) :: paths
    type(mesh_options) :: mesh
    type(io_options) :: io
    type(variable_options) :: variables
    type(solver_options) :: solve

    !! Parallel-related
    logical :: use_mpi_splitting
    integer :: split_type
  end type ccs_options


  interface
    !v Subroutine to get the runtime configuration.
    module subroutine get_config(par_env, run_options)
      class(parallel_environment), intent(in) :: par_env
      type(ccs_options), intent(out) :: run_options
    end subroutine get_config

    !v Subroutine to configure sub parallel environments.
    !
    ! Currently this only configures the shared environment, but it could be extended for particle
    ! environments, etc.
    module subroutine configure_parallelism(run_options, par_env, shared_env)
      type(ccs_options), intent(in) :: run_options
      class(parallel_environment), allocatable, intent(in) ::  par_env
      class(parallel_environment), allocatable, intent(out) :: shared_env
    end subroutine configure_parallelism
    

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
