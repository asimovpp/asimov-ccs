!v Module file core.mod
!
! The core module implements the core functionality to define and run a CCS case.

module core

  use kinds, only: ccs_int, ccs_real
  use types, only: fluid

  use parallel_types, only: parallel_environment
  use constants, only: ccs_string_len


  implicit none

  private

  public :: initialise_flow
  public :: configure_parallelism
  public :: initialise_mesh
  public :: initialise_fields
  public :: run_solver
  public :: ccs_options
  integer(ccs_int), parameter, public :: read_input_mesh = 1
  integer(ccs_int), parameter, public :: build_mesh_2d = 2
  integer(ccs_int), parameter, public :: build_mesh_3d = 3

  !v Type to contain the configuration of the CCS run
  type :: ccs_options
    !! Generic/misc
    character(len=:), allocatable :: case_name
    character(len=:), allocatable :: case_path
    character(len=:), allocatable :: config_file
    character(len = ccs_string_len), dimension(:), allocatable :: variable_names  
    integer(ccs_int), dimension(:), allocatable :: variable_types              

    !! Parallel-related
    logical :: use_mpi_splitting
    integer :: split_type

    !! Mesh-related
    character(len=:), allocatable :: mesh_path

    !! Solver-related
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

    module subroutine initialise_fields(par_env, run_options, flow_fields)
      class(parallel_environment), intent(in), allocatable:: par_env    !< The parallel environment
      type(ccs_options), intent(in) :: run_options    !< Object containing relevant options for building/reading the mesh
      type(fluid), intent(out) :: flow_fields         !< The fluid fields object being initialised
    end subroutine initialise_fields

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
