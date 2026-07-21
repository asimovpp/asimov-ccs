!v Module file core.mod
!
! The core module implements the core functionality to define and run a CCS case.

module core

  use constants, only: ccs_string_len
  use kinds, only: ccs_int, ccs_real
  use types, only: fluid, solver_params

  use parallel_types, only: parallel_environment
  use constants, only: ccs_string_len

  implicit none

  private

  public :: get_config
  public :: initialise_flow
  public :: configure_parallelism
  public :: initialise_mesh
  public :: initialise_fields
  public :: run_solver
  public :: ccs_options

  integer(ccs_int), parameter :: mesh_null = 0
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
    integer(ccs_int) :: init_mesh_type = mesh_null
    integer(ccs_int) :: cps = huge(0_ccs_int)
    real(ccs_real) :: domain_size = huge(0.0_ccs_real)
    logical :: compute_bwidth = .true.
    logical :: compute_partqual = .true.
    character(len=ccs_string_len), dimension(:), allocatable :: bnd_names
  end type mesh_options
  
  !v Options for IO configuration
  type :: io_options
    integer(ccs_int) :: write_frequency = huge(0_ccs_int)
    logical :: write_gradients = .false.
  end type io_options

  !v Options for variable declarations
  type :: variable_options
    character(len=ccs_string_len), dimension(:), allocatable :: variable_names
    integer(ccs_int), dimension(:), allocatable :: variable_types
    character(len=ccs_string_len), dimension(:), allocatable :: output_variables
    logical :: restart = .false.
  end type variable_options

  !v Options for solver configuration
  type :: solver_options
    logical :: unsteady = .false.
    integer(ccs_int) :: num_steps = huge(0_ccs_int)
    integer(ccs_int) :: num_iters = huge(0_ccs_int)
    integer(ccs_int) :: it_start = huge(0_ccs_int)
    integer(ccs_int) :: it_end = huge(0_ccs_int)
    real(ccs_real) :: dt = huge(0.0_ccs_real)
    type(solver_params), dimension(:), allocatable :: solver_eq_parameters
    real(ccs_real) :: global_res_target = huge(0.0_ccs_real)
    logical :: debug = .false. ! Turn on solver checks
  end type solver_options

  !v Options for parallelism
  type :: parallel_options
    logical :: use_mpi_splitting = .false.
    integer :: split_type = huge(0_ccs_int)
  end type parallel_options
  
  !v Reference values for the problem
  type :: ref_vals
    real(ccs_real) :: p_ref = huge(0.0_ccs_real)    !< reference pressure
    real(ccs_real) :: p_total = huge(0.0_ccs_real)  !< total pressure
    real(ccs_real) :: temp_ref = huge(0.0_ccs_real) !< reference temperature
    real(ccs_real) :: dens_ref = huge(0.0_ccs_real) !< reference density
    real(ccs_real) :: visc_ref = huge(0.0_ccs_real) !< laminar viscosity
    real(ccs_real) :: velo_ref = huge(0.0_ccs_real) !< reference velocity
    real(ccs_real) :: len_ref = huge(0.0_ccs_real)  !< reference length, used to define the Reynolds number of the flow
    integer :: pref_at_cell = huge(0_ccs_int)       !< cell at which the reference pressure is set
  end type ref_vals
  
  !v Type to contain the configuration of the CCS run
  type :: ccs_options
    type(ccs_paths) :: paths
    type(mesh_options) :: mesh
    type(io_options) :: io
    type(variable_options) :: variables
    type(solver_options) :: solve
    type(parallel_options) :: parallel
    type(ref_vals) :: reference_values
  end type ccs_options

  interface
    !v Subroutine to get the runtime configuration.
    module subroutine get_config(par_env, run_options)
      class(parallel_environment), intent(in) :: par_env !< The parallel environment
      type(ccs_options), intent(out) :: run_options      !< The runtime configuration
    end subroutine get_config

    !v Subroutine to configure sub parallel environments.
    !
    ! Currently this only configures the shared environment, but it could be extended for particle
    ! environments, etc.
    module subroutine configure_parallelism(run_options, par_env, shared_env)
      type(ccs_options), intent(in) :: run_options
      class(parallel_environment), allocatable, intent(in) ::  par_env    !< The (primary) parallel environment
      class(parallel_environment), allocatable, intent(out) :: shared_env !< The split (shared) parallel environment
    end subroutine configure_parallelism
    
    !v Initialise both cell centre values and mass fluxes by calling get_init_flow and get_init_mass_flux
    !  on every cell or face
    module subroutine initialise_flow(par_env, run_options, flow_fields, get_init_flow, get_init_mass_flux)
      class(parallel_environment), intent(in) :: par_env !< Parallel environment
      type(ccs_options), intent(in) :: run_options       !< Runtime configuration
      type(fluid), intent(inout) :: flow_fields          !< The flow
      interface 
        !> User-supplied subroutine to set field values at cell centres
        pure subroutine get_init_flow(loc_p, field_name, init_val)
          use kinds, only: ccs_real
          use types, only: cell_locator
          type(cell_locator), intent(in) :: loc_p    !< Cell locator
          character(len=*), intent(in) :: field_name !< Name of field being initialised
          real(ccs_real), intent(inout) :: init_val  !< The initial value to set for the field
        end subroutine

        !> User-supplied subroutine to set the mass flux at face centres
        pure subroutine get_init_mass_flux(loc_f, init_val)
          use kinds, only: ccs_real
          use types, only: face_locator
          type(face_locator), intent(in) :: loc_f   !< Face locator
          real(ccs_real), intent(inout) :: init_val !< The initial value to set for the mass flux
        end subroutine
      end interface
    end subroutine

    !v Subroutine to initialise the mesh.
    !
    ! This is responsible for the building or reading of the mesh (whichever is specified in the config file)
    module subroutine initialise_mesh(par_env, shared_env, run_options)
      class(parallel_environment), intent(in), allocatable :: par_env    !< The (primary) parallel environment
      class(parallel_environment), intent(in), allocatable :: shared_env !< The split (shared) parallel environment
      type(ccs_options), intent(in) :: run_options
    end subroutine initialise_mesh

    !v Subroutine to initialise fields
    !
    ! This is responsible for setting the building common fields and any other fields specified in the case config file
    module subroutine initialise_fields(par_env, run_options, flow_fields)
      class(parallel_environment), intent(in), allocatable:: par_env !< The parallel environment
      type(ccs_options), intent(in) :: run_options                   !< The runtime configuration
      type(fluid), intent(out) :: flow_fields                        !< The flow field structure
    end subroutine initialise_fields

    !v Subroutine to run a flow problem.
    !
    ! This is responsible for the time loop, calling post-processing subroutines and performing
    ! solution output.
    module subroutine run_solver(par_env, run_options, eval_sources, postproc, flow_fields)
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
        end subroutine postproc
      end interface
      type(fluid), intent(inout) :: flow_fields !< The flow field structure, contains the solution
    end subroutine run_solver

  end interface
    
end module core
