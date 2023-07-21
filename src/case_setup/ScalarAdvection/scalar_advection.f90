!v Program file for scalar advection case

program scalar_advection

  ! ASiMoV-CCS uses
  use kinds, only: ccs_real, ccs_int
  use case_config, only: velocity_solver_method_name, velocity_solver_precon_name, &
                         pressure_solver_method_name, pressure_solver_precon_name
  use types, only: vector_spec, ccs_vector, matrix_spec, ccs_matrix, &
                   equation_system, linear_solver, ccs_mesh, &
                   field, upwind_field, central_field, bc_config
  use constants, only: ccs_split_type_shared, ccs_split_type_low_high, ccs_split_undefined
  use vec, only: create_vector
  use mat, only: create_matrix, set_nnz
  use solver, only: create_solver, solve, set_equation_system
  use utils, only: update, initialise, set_size
  use mesh_utils, only: build_square_mesh
  use parallel_types, only: parallel_environment
  use parallel, only: initialise_parallel_environment, create_new_par_env &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments
  use fv, only: compute_fluxes
  use boundary_conditions, only: read_bc_config

  implicit none

  class(parallel_environment), allocatable, target :: par_env
  class(parallel_environment), allocatable, target :: shared_env
  class(ccs_vector), allocatable, target :: source
  class(ccs_vector), allocatable :: solution
  class(ccs_matrix), allocatable, target :: M
  class(linear_solver), allocatable :: scalar_solver

  type(vector_spec) :: vec_properties
  type(matrix_spec) :: mat_properties
  type(equation_system) :: scalar_equation_system
  type(ccs_mesh) :: mesh

  class(field), allocatable :: mf !< Prescribed face velocity field
  class(field), allocatable :: scalar

  integer(ccs_int) :: cps = 50 ! Default value for cells per side
  integer(ccs_int) :: direction = 0 ! pass zero for "direction" of scalar field when computing fluxes

  double precision :: start_time
  double precision :: end_time
  logical :: use_mpi_splitting

  call initialise_parallel_environment(par_env)
  use_mpi_splitting = .false.
  call create_new_par_env(par_env, ccs_split_type_low_high, use_mpi_splitting, shared_env)

  ! set solver and preconditioner info
  velocity_solver_method_name = "gmres"
  velocity_solver_precon_name = "bjacobi"
  pressure_solver_method_name = "cg"
  pressure_solver_precon_name = "gamg"

  call read_command_line_arguments(par_env)
  call timer(start_time)

  ! Set up the square mesh
  mesh = build_square_mesh(par_env, shared_env, cps, 1.0_ccs_real)

  ! Init velocities and scalar
  allocate (central_field :: mf)
  allocate (upwind_field :: scalar)

  ! Read bc configuration
  call read_bc_config("./case_setup/ScalarAdvection/ScalarAdvection_config.yaml", "scalar", scalar)

  ! Initialise with default values
  call initialise(vec_properties)
  call initialise(mat_properties)
  call initialise(scalar_equation_system)

  ! Create stiffness matrix
  call set_size(par_env, mesh, mat_properties)
  call set_nnz(5, mat_properties)
  call create_matrix(mat_properties, M)

  ! Create right-hand-side and solution vectors
  call set_size(par_env, mesh, vec_properties)
  call create_vector(vec_properties, source)
  call create_vector(vec_properties, solution)
  call create_vector(vec_properties, scalar%values)
  call create_vector(vec_properties, mf%values)

  ! Set advection velocity
  call set_advection_velocity(mesh, mf)

  ! Actually compute the values to fill the matrix
  call compute_fluxes(scalar, mf, mesh, direction, M, source)

  call update(M) ! parallel assembly for M

  call update(source) ! parallel assembly for source

  ! Create linear solver & set options
  call set_equation_system(par_env, source, scalar%values, M, scalar_equation_system)
  call create_solver(scalar_equation_system, scalar_solver)
  call solve(scalar_solver)

  ! Clean up
  deallocate (scalar)
  deallocate (source)
  deallocate (solution)
  deallocate (M)
  deallocate (mf)
  deallocate (scalar_solver)

  call timer(end_time)

  if (par_env%proc_id == par_env%root) then
    print *, "Elapsed time = ", (end_time - start_time)
  end if

  call cleanup_parallel_environment(par_env)

contains

  subroutine set_advection_velocity(mesh, mf)
    use constants, only: add_mode
    use types, only: vector_values, cell_locator
    use meshing, only: create_cell_locator, get_global_index, get_local_num_cells
    use fv, only: calc_cell_coords
    use utils, only: set_row, set_entry, set_values

    class(ccs_mesh), intent(in) :: mesh
    class(field), intent(inout) :: mf
    integer(ccs_int) :: row, col
    integer(ccs_int) :: index_p, global_index_p, n_local
    real(ccs_real) :: mf_val
    type(cell_locator) :: loc_p
    type(vector_values) :: mf_vals

    real(ccs_real) :: u, v

    mf_vals%setter_mode = add_mode

    call get_local_num_cells(mesh, n_local)

    allocate (mf_vals%global_indices(n_local))
    allocate (mf_vals%values(n_local))

    ! Set IC velocity and scalar fields
    do index_p = 1, n_local
      call create_cell_locator(mesh, index_p, loc_p)
      call get_global_index(loc_p, global_index_p)
      call calc_cell_coords(global_index_p, cps, row, col)

      ! TODO: this should be in a face loop, compute these based on normals and set mf appropriately
      u = real(col, ccs_real) / real(cps, ccs_real)
      v = -real(row, ccs_real) / real(cps, ccs_real)

      mf_val = u + v

      call set_row(global_index_p, mf_vals)
      call set_entry(mf_val, mf_vals)
    end do
    call set_values(mf_vals, mf%values)

    deallocate (mf_vals%global_indices)
    deallocate (mf_vals%values)

  end subroutine set_advection_velocity

end program scalar_advection
