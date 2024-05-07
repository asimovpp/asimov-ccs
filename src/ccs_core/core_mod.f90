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
  public :: initialise_flow
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
    integer(ccs_int) :: cps = huge(0)
    real(ccs_real) :: domain_size
    logical :: compute_bwidth = .true.
    logical :: compute_partqual = .true.
  end type mesh_options
  
  !v Options for IO configuration
  type :: io_options
    integer(ccs_int) :: write_frequency = huge(0)
  end type io_options

  !v Options for variable declarations
  type :: variable_options
    character(len=ccs_string_len), dimension(:), allocatable :: variable_names
    integer(ccs_int), dimension(:), allocatable :: variable_types
  end type variable_options

  !v Options for solver configuration
  type :: solver_options
    integer(ccs_int) :: num_steps = huge(0)
    integer(ccs_int) :: num_iters = huge(0)
    integer(ccs_int) :: it_start
    integer(ccs_int) :: it_end
    real(ccs_real) :: dt = huge(0.0_ccs_real)
    real(ccs_real) :: res_target = huge(0.0_ccs_real)
    real(ccs_real) :: velocity_relax = huge(0.0_ccs_real)
    real(ccs_real) :: pressure_relax = huge(0.0_ccs_real)
  end type solver_options

  !v Options for parallelism
  type :: parallel_options
    logical :: use_mpi_splitting
    integer :: split_type
  end type parallel_options
  
  !v Type to contain the configuration of the CCS run
  type :: ccs_options
    type(ccs_paths) :: paths
    type(mesh_options) :: mesh
    type(io_options) :: io
    type(variable_options) :: variables
    type(solver_options) :: solve
    type(parallel_options) :: parallel
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
    
    !> Initialise both cell centre values and mass fluxes by calling get_init_flow and get_init_mass_flux
    !  on every cell or face
    module subroutine initialise_flow(flow_fields, get_init_flow, get_init_mass_flux)
      type(fluid), intent(inout) :: flow_fields
      interface 
        pure subroutine get_init_flow(loc_p, field_name, init_val)
          use kinds, only: ccs_real
          use types, only: cell_locator
          type(cell_locator), intent(in) :: loc_p
          character(len=*), intent(in) :: field_name
          real(ccs_real), intent(inout) :: init_val
        end subroutine

        pure subroutine get_init_mass_flux(loc_f, init_val)
          use kinds, only: ccs_real
          use types, only: face_locator
          type(face_locator), intent(in) :: loc_f
          real(ccs_real), intent(inout) :: init_val
        end subroutine
      end interface
    end subroutine

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
    module subroutine run_solver(par_env, run_options, postproc, flow_fields)
      class(parallel_environment), allocatable, intent(in) :: par_env
      type(ccs_options), intent(in) :: run_options
      interface
        subroutine postproc(par_env, flow_fields)
          use types, only: fluid
          use parallel_types, only: parallel_environment
          class(parallel_environment), allocatable, intent(in) :: par_env
          type(fluid), intent(in) :: flow_fields
        end subroutine postproc
      end interface
      type(fluid), intent(inout) :: flow_fields
    end subroutine run_solver

  end interface
    
end module core
