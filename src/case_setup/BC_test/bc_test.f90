!>  Program file for scalar advection case
!
!

program scalar_advection

  !! ASiMoV-CCS uses
  use kinds, only : ccs_real, ccs_int
  use types, only : vector_spec, ccs_vector, matrix_spec, ccs_matrix, &
                    equation_system, linear_solver, ccs_mesh, &
                    field, upwind_field, central_field, bc_config
  use vec, only : create_vector
  use mat, only : create_matrix, set_nnz
  use solver, only : create_solver, solve, set_equation_system
  use utils, only : update, initialise, set_size
  use mesh_utils, only : build_square_mesh
  use parallel_types, only: parallel_environment
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments
  use fv, only : compute_fluxes
  use boundary_conditions, only : read_bc_config

  implicit none

  class(parallel_environment), allocatable, target :: par_env
  class(ccs_vector), allocatable, target :: source
  class(ccs_vector), allocatable :: solution
  class(ccs_matrix), allocatable, target :: M
  class(linear_solver), allocatable :: scalar_solver

  type(vector_spec) :: vec_properties
  type(matrix_spec) :: mat_properties
  type(equation_system) :: scalar_equation_system
  type(ccs_mesh) :: mesh

  class(field), allocatable :: mf          ! Prescribed face velocity field
  class(field), allocatable :: scalar

  integer(ccs_int) :: cps = 50 ! Default value for cells per side

  double precision :: start_time
  double precision :: end_time

  call initialise_parallel_environment(par_env) 
  call read_command_line_arguments(par_env)
  call timer(start_time)

  ! Set up the square mesh
  mesh = build_square_mesh(par_env, cps, 1.0_ccs_real)

  ! Init velocities and scalar
  allocate(central_field :: mf)
  allocate(upwind_field :: scalar)

  ! Read bc configuration
  call read_bc_config("BC_test_config.yaml", mf%bcs)

  call cleanup_parallel_environment(par_env)

end program scalar_advection
