!> @brief Program file for scalar advection case
!
!

program scalar_advection

  !! ASiMoV-CCS uses
  use kinds, only : accs_real, accs_int
  use types, only : vector_init_data, vector, matrix_init_data, ccs_matrix, &
                    linear_system, linear_solver, mesh, &
                    field, upwind_field, central_field, bc_config
  use vec, only : create_vector
  use mat, only : create_matrix, set_nnz
  use solver, only : create_solver, solve, set_linear_system
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
  class(vector), allocatable, target :: source
  class(vector), allocatable :: solution
  class(ccs_matrix), allocatable, target :: M
  class(linear_solver), allocatable :: scalar_solver

  type(vector_init_data) :: vec_sizes
  type(matrix_init_data) :: mat_sizes
  type(linear_system) :: scalar_linear_system
  type(mesh) :: square_mesh
  type(bc_config) :: bcs  !XXX: BCs are part of the fields structure now. fix this.

  class(field), allocatable :: mf          ! Prescribed face velocity field
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
  square_mesh = build_square_mesh(par_env, cps, 1.0_accs_real)

  ! Init velocities and scalar
  allocate(central_field :: mf)
  allocate(upwind_field :: scalar)

  ! Initialise with default values
  call initialise(mat_sizes)
  call initialise(vec_sizes)
  call initialise(scalar_linear_system)

  ! Create stiffness matrix
  call set_size(par_env, square_mesh, mat_sizes)
  call set_nnz(5, mat_sizes) 
  call create_matrix(mat_sizes, M)

  ! Create right-hand-side and solution vectors
  call set_size(par_env, square_mesh, vec_sizes)
  call create_vector(vec_sizes, source)
  call create_vector(vec_sizes, solution)
  call create_vector(vec_sizes, scalar%vec)
  call create_vector(vec_sizes, mf%vec)

  ! Set advection velocity
  call set_advection_velocity(square_mesh, mf)

  ! Actually compute the values to fill the matrix
  call compute_fluxes(scalar, mf, square_mesh, bcs, cps, M, source)

  call update(M) ! parallel assembly for M

  call update(source) ! parallel assembly for source

  ! Create linear solver & set options
  call set_linear_system(par_env, source, scalar%vec, M, scalar_linear_system)
  call create_solver(scalar_linear_system, scalar_solver)
  call solve(scalar_solver)
  
  ! Clean up
  deallocate(scalar)
  deallocate(source)
  deallocate(solution)
  deallocate(M)
  deallocate(mf)
  deallocate(scalar_solver)

  call timer(end_time)

  if (par_env%proc_id == par_env%root) then
    print *, "Elapsed time = ", (end_time - start_time)
  end if

  call cleanup_parallel_environment(par_env)

contains

  subroutine set_advection_velocity(cell_mesh, mf)
    use constants, only: add_mode
    use types, only: vector_values, cell_locator
    use meshing, only: set_cell_location, get_global_index
    use fv, only: calc_cell_coords
    use utils, only: pack_entries, set_values

    class(mesh), intent(in) :: cell_mesh
    class(field), intent(inout) :: mf
    integer(accs_int) :: row, col
    integer(accs_int) :: local_idx, self_idx
    real(accs_real) :: mf_val
    type(cell_locator) :: self_loc
    type(vector_values) :: mf_vals

    real(accs_real) :: u, v

    mf_vals%mode = add_mode

    associate(n_local => cell_mesh%nlocal)
      allocate(mf_vals%idx(n_local))
      allocate(mf_vals%val(n_local))
      
      ! Set IC velocity and scalar fields
      do local_idx = 1, n_local
        call set_cell_location(cell_mesh, local_idx, self_loc)
        call get_global_index(self_loc, self_idx)
        call calc_cell_coords(self_idx, cps, row, col)

        ! TODO: this should be in a face loop, compute these based on normals and set mf appropriately
        u = real(col, accs_real)/real(cps, accs_real)
        v = -real(row, accs_real)/real(cps, accs_real)

        mf_val = u + v
        
        call pack_entries(local_idx, self_idx, mf_val, mf_vals)
      end do
    end associate
    call set_values(mf_vals, mf%vec)

    deallocate(mf_vals%idx, mf_vals%val)
  end subroutine set_advection_velocity

end program scalar_advection
