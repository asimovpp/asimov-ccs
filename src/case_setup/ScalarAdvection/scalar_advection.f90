!> @brief Program file for scalar advection case
!
!

program scalar_advection

  !! ASiMoV-CCS uses
  use kinds, only : accs_real, accs_int
  use types, only : vector_init_data, vector, matrix_init_data, matrix, &
                    linear_system, linear_solver, mesh, set_global_matrix_size, &
                    field, upwind_field, central_field, bc_config
  !use vec, only : create_vector
  use vec, only : create_vector, vec_view
  use mat, only : create_matrix, set_nnz
  use solver, only : create_solver, solve, set_linear_system, axpy, norm
  use utils, only : update, begin_update, end_update, finalise, initialise, &
                    set_global_size
  use mesh_utils, only : build_square_mesh
  use petsctypes, only : matrix_petsc
  use parallel_types, only: parallel_environment
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments
  use fv, only : compute_fluxes
  use boundary_conditions, only : read_bc_config

  implicit none

  class(parallel_environment), allocatable, target :: par_env
  class(vector), allocatable, target :: source
  class(vector), allocatable :: solution
  class(matrix), allocatable, target :: M
  class(linear_solver), allocatable :: scalar_solver

  type(vector_init_data) :: vec_sizes
  type(matrix_init_data) :: mat_sizes
  type(linear_system) :: scalar_linear_system
  type(mesh) :: square_mesh
  type(bc_config) :: bcs

  class(field), allocatable :: u, v          ! Prescribed x, y velocity fields
  class(field), allocatable :: scalar

  integer(accs_int) :: cps = 50 ! Default value for cells per side

  double precision :: start_time
  double precision :: end_time

  call initialise_parallel_environment(par_env) 
  call read_command_line_arguments(par_env)
  call timer(start_time)

  ! Read bc configuration
  call read_bc_config("./case_setup/ScalarAdvection/ScalarAdvection_config.yaml", bcs)

  ! Set up the square mesh
  square_mesh = build_square_mesh(cps, 1.0_accs_real, par_env)

  ! Init velocities and scalar
  allocate(upwind_field :: u)
  allocate(upwind_field :: v)
  allocate(upwind_field :: scalar)

  ! Initialise with default values
  call initialise(mat_sizes)
  call initialise(vec_sizes)
  call initialise(scalar_linear_system)

  ! Create stiffness matrix
  call set_global_size(mat_sizes, square_mesh%nglobal, square_mesh%nglobal, par_env)
  call set_nnz(mat_sizes, 5) 
  call create_matrix(mat_sizes, M)

  ! Create right-hand-side and solution vectors
  call set_global_size(vec_sizes, square_mesh%nglobal, par_env)
  call create_vector(vec_sizes, source)
  call create_vector(vec_sizes, solution)
  call create_vector(vec_sizes, scalar%vec)
  call create_vector(vec_sizes, u%vec)
  call create_vector(vec_sizes, v%vec)

  ! Set advection velocity
  call set_advection_velocity(square_mesh, u, v)

  ! Actually compute the values to fill the matrix
  call compute_fluxes(u, v, square_mesh, bcs, cps, M, source)

  call update(M) ! parallel assembly for M

  call update(source) ! parallel assembly for source

  ! Create linear solver & set options
  call set_linear_system(scalar_linear_system, source, scalar%vec, M, par_env)
  call create_solver(scalar_linear_system, scalar_solver)
  call solve(scalar_solver)

  call vec_view(scalar%vec)
  
  ! Clean up
  deallocate(scalar)
  deallocate(source)
  deallocate(solution)
  deallocate(M)
  deallocate(u)
  deallocate(v)
  deallocate(scalar_solver)

  call timer(end_time)

  if (par_env%proc_id == par_env%root) then
    print *, "Elapsed time = ", (end_time - start_time)
  end if

  call cleanup_parallel_environment(par_env)

contains

  subroutine set_advection_velocity(cell_mesh, u, v)
    use constants, only: add_mode
    use types, only: vector_values, cell_locator
    use meshing, only: set_cell_location, get_global_index
    use fv, only: calc_cell_coords
    use utils, only: pack_entries, set_values
    class(mesh), intent(in) :: cell_mesh
    class(field), intent(inout) :: u, v
    integer(accs_int) :: i
    integer(accs_int) :: row, col
    integer(accs_int) :: local_idx, self_idx
    real(accs_real) :: u_val, v_val
    type(cell_locator) :: self_loc
    type(vector_values) :: u_vals, v_vals

    u_vals%mode = add_mode
    v_vals%mode = add_mode

    associate(n_local => cell_mesh%nlocal)
      allocate(u_vals%idx(n_local))
      allocate(v_vals%idx(n_local))
      allocate(u_vals%val(n_local))
      allocate(v_vals%val(n_local))
      
      ! Set IC velocity and scalar fields
      do local_idx = 1, n_local
        call set_cell_location(self_loc, cell_mesh, local_idx)
        call get_global_index(self_loc, self_idx)
        call calc_cell_coords(self_idx, cps, row, col)

        u_val = real(col, accs_real)/real(cps, accs_real)
        v_val = -real(row, accs_real)/real(cps, accs_real)

        call pack_entries(u_vals, local_idx, self_idx, u_val)
        call pack_entries(v_vals, local_idx, self_idx, v_val)
      end do
    end associate
    call set_values(u_vals, u%vec)
    call set_values(v_vals, v%vec)

    deallocate(u_vals%idx, v_vals%idx, u_vals%val, v_vals%val)
  end subroutine set_advection_velocity

end program scalar_advection
