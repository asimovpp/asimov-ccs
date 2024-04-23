!v Module file core.mod
!
! The core module implements the core functionality to define and run a CCS case.

module core

  use kinds, only: ccs_int, ccs_real
  use types, only: fluid

  use parallel_types, only: parallel_environment


  implicit none

  private

  public :: core_initialise_flow
  public :: core_initialise_mass_flux
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

    !> Initialise every cell based field prompting values from get_init_flow
    module subroutine core_initialise_flow(flow_fields, get_init_flow)
      type(fluid), intent(inout) :: flow_fields
      interface 
        pure subroutine get_init_flow(loc_p, field_name, init_val)
          use kinds, only: ccs_real
          use types, only: cell_locator
          type(cell_locator), intent(in) :: loc_p
          character(len=*), intent(in) :: field_name
          real(ccs_real), intent(inout) :: init_val
        end subroutine
      end interface
    end subroutine

    !> Initialise mass flux field using calling get_init_mass_flux
    module subroutine core_initialise_mass_flux(flow_fields, get_init_mass_flux)
      type(fluid), intent(inout) :: flow_fields
      interface
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
    module subroutine run_solver(par_env, run_options, flow_fields)
      class(parallel_environment), allocatable, intent(in) :: par_env
      type(ccs_options), intent(in) :: run_options
      type(fluid), intent(inout) :: flow_fields
    end subroutine run_solver

  end interface
  
end module core
