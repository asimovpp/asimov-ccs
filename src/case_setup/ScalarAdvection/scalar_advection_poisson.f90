!> @brief Program file for scalar advection case
!
!

program scalar_advection

  !! ASiMoV-CCS uses
  use kinds, only : accs_real, accs_int
  use types, only : vector_init_data, vector, matrix_init_data, matrix, &
                    linear_system, linear_solver, mesh, set_global_matrix_size, &
                    field, upwind_field, central_field
  use vec, only : create_vector, axpy, norm
  use mat, only : create_matrix, set_nnz
  use solver, only : create_solver, solve, set_linear_system
  use utils, only : update, begin_update, end_update, finalise, initialise, &
                    set_global_size
  use mesh_utils, only : build_square_mesh
  use petsctypes, only : matrix_petsc
  use parallel_types, only: parallel_environment
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments
  use fv, only : compute_fluxes

  use petsc, only: ADD_VALUES  

  implicit none

  class(parallel_environment), allocatable, target :: par_env
  class(vector), allocatable, target :: scalar, source
  class(vector), allocatable :: solution
  class(matrix), allocatable, target :: M
  class(linear_solver), allocatable :: scalar_solver

  type(vector_init_data) :: vec_sizes
  type(matrix_init_data) :: mat_sizes
  type(linear_system) :: scalar_linear_system
  type(mesh) :: square_mesh

  class(field), allocatable :: u, v          ! Prescribed x, y velocity fields

  integer(accs_int) :: cps = 50 ! Default value for cells per side

  double precision :: start_time
  double precision :: end_time

  call initialise_parallel_environment(par_env) 
  call read_command_line_arguments(par_env)
  call timer(start_time)

  ! Init ICs (velocities, BC scalar, mesh, etc)
  allocate(upwind_field :: u)
  allocate(upwind_field :: v)
  call initialise_scalar_advection(par_env, u, v)

  !! Initialise with default values
  call initialise(mat_sizes)
  call initialise(vec_sizes)
  call initialise(scalar_linear_system)

  !! Create stiffness matrix
  call set_global_size(mat_sizes, square_mesh%n, square_mesh%n, par_env)
  call set_nnz(mat_sizes, 5) 
  call create_matrix(mat_sizes, M)

  !! Create right-hand-side and solution vectors
  call set_global_size(vec_sizes, square_mesh%n, par_env)
  call create_vector(vec_sizes, source)
  call create_vector(vec_sizes, solution)
  call create_vector(vec_sizes, scalar)

  ! Actually compute the values to fill the matrix
  call compute_fluxes(u, v, square_mesh, cps, M, source)

  call begin_update(M) ! Start the parallel assembly for M

  call begin_update(source) ! Start the parallel assembly for source
  call end_update(M) ! Complete the parallel assembly for M
  call end_update(source) ! Complete the parallel assembly for source

  !! Create linear solver & set options
  call set_linear_system(scalar_linear_system, source, scalar, M, par_env)
  call create_solver(scalar_linear_system, scalar_solver)
  call solve(scalar_solver)
  
  !! Clean up
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

  ! Make the square mesh, and setup ICs for u, v, and our scalar field
  subroutine initialise_scalar_advection(par_env, u, v)

    class(parallel_environment) :: par_env
    class(field), intent(inout) :: u, v

    integer(accs_int) :: i

    square_mesh = build_square_mesh(cps, 1.0_accs_real, par_env)

    ! Allocate velocity and scalar field arrays
    allocate(u%val(cps,cps))
    allocate(v%val(cps,cps))

    ! Set IC velocity and scalar fields
    do i = 1, cps
      u%val(i,:) = real(i, accs_real)/real(cps, accs_real)
      v%val(:,i) = -real(i, accs_real)/real(cps, accs_real)
    end do
  end subroutine initialise_scalar_advection

end program scalar_advection
