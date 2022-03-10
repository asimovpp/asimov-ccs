!> @brief Program file for pressure-velocity coupling case
!
!> @details This case demonstrates solution of the Navier-Stokes equations
!!          using the SIMPLE algorithm for pressure-velocity coupling.
!

program simple

  use kinds, only: accs_real, accs_int
  use types, only: field, upwind_field, central_field, mesh, &
                   vector_init_data, matrix, vector
  use parallel, only: initialise_parallel_environment
  use parallel_types, only: parallel_environment
  use mesh_utils, only: build_square_mesh
  use vec, only: create_vector
  use petsctypes, only: matrix_petsc, vector_petsc
  use pv_coupling, only: solve_nonlinear
  use mat, only: get_matrix_diagonal
                      
  implicit none

  class(parallel_environment), allocatable, target :: par_env
  type(mesh)             :: square_mesh
  type(vector_init_data) :: vec_sizes

  class(field), allocatable :: u, v, p, pp

  integer(accs_int) :: cps = 50 ! Default value for cells per side

  integer(accs_int) :: it_start, it_end

  double precision :: start_time
  double precision :: end_time
 
  ! Set start and end iteration numbers (eventually will be read from input file)
  it_start = 1
  it_end   = 10

  call initialise_parallel_environment(par_env)
  call read_command_line_arguments(par_env)

  call sync(par_env)
  call timer(start_time)

  ! Create a square mesh
  square_mesh = build_square_mesh(cps, 1.0_accs_real, par_env)

  ! Initialise fields
  allocate(upwind_field :: u)
  allocate(upwind_field :: v)
  allocate(central_field :: p)
  allocate(central_field :: pp)

  ! Create and initialise field vectors
  call set_global_size(vec_sizes, square_mesh%nglobal, par_env)
  call create_vector(vec_sizes, u%vec)
  call create_vector(vec_sizes, v%vec)
  call create_vector(vec_sizes, p%vec)
  call create_vector(vec_sizes, pp%vec)

  ! Initialise velocity field
  call initialise_velocity(square_mesh, u, v)

  ! Solve using SIMPLE algorithm
  call solve_nonlinear(par_env, square_mesh, it_start, it_end, u, v, p, pp)

  ! Clean-up
  deallocate(u)
  deallocate(v)
  deallocate(p)
  deallocate(pp)

  call timer(end_time)

  if (par_env%proc_id == par_env%root) then
    print *, "Elapsed time = ", (end_time - start_time)
  end if

  call cleanup_parallel_environment(par_env)

contains

  subroutine initialise_velocity(cell_mesh, u, v)

    use constants, only: add_mode
    use types, only: vector_values, cell_locator
    use meshing, only: set_cell_location, get_global_index
    use fv, only: calc_cell_coords
    use utils, only: pack_entries, set_values

    ! Arguments
    class(mesh), intent(in) :: cell_mesh
    class(field), intent(inout) :: u, v

    ! Local variables
    integer(accs_int) :: i
    integer(accs_int) :: row, col
    integer(accs_int) :: local_idx, self_idx
    real(accs_real) :: u_val, v_val
    type(cell_locator) :: self_loc
    type(vector_values) :: u_vals, v_vals

    ! Set mode
    u_vals%mode = add_mode
    v_vals%mode = add_mode

    ! Set alias
    associate(n_local => cell_mesh%nlocal)
      ! Allocate temporary arrays for storing global cell indices 
      allocate(u_vals%idx(n_local))
      allocate(v_vals%idx(n_local))

      ! Allocate temporary arrays for storing values
      allocate(u_vals%val(n_local))
      allocate(v_vals%val(n_local))

      ! Set initial values for velocity fields
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

  end subroutine initialise_velocity

end program simple