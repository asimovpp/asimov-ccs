!> @brief Program file for pressure-velocity coupling case
!
!> @details This case demonstrates solution of the Navier-Stokes equations
!!          using the SIMPLE algorithm for pressure-velocity coupling.
!

program simple

  use kinds, only: accs_real, accs_int
  use parallel, only: initialise_parallel_environment, &
                      


  implicit none

  class(parallel_environment), allocatable, target :: par_env
  class(vector), allocatable, target :: source
  class(vector), allocatable :: solution
  class(matrix), allocatable, target :: M
  class(linear_solver), allocatable :: scalar_solver

  type(vector_init_data) :: vec_sizes
  type(matrix_init_data) :: mat_sizes
  type(linear_system)    :: scalar_linear_system
  type(mesh)             :: square_mesh
  type(bc_config)        :: bcs

  class(field), allocatable :: u, v, p

  integer(accs_int) :: cps = 50 ! Default value for cells per side

  double precision :: start_time
  double precision :: end_time
    
  call initialise_parallel_environment(par_env)
  call read_command_line_arguments(par_env)
  call timer(start_time)

  ! Initialise fields
  allocate(upwind_field :: u)
  allocate(upwind_field :: v)
  allocate(central_field :: p)
  call initialise_fields(par_env, u, v, p)

  ! Initialise linear system
  call initialise(mat_sizes)
  call initialise(vec_sizes)
  call initialise(scalar_linear_system)

  ! Create coefficient matrix
  call set_global_size(mat_sizes, square_mesh%n, square_mesh%n, par_env)
  call set_nnz(mat_sizes, 5)
  call create_matrix(mat_sizes, M)

  ! Create RHS and solution vectors
  call set_global_size(vec_sizes, square_mesh%n, par_env)
  call create_vector(vec_sizes, source)
  call create_vector(vec_sizes, solution)

  ! Solve using SIMPLE algorithm
  call solve_nonlinear(u, v, p, M, solution, source, square_mesh)

  ! Clean-up
  deallocate(source)
  deallocate(solution)
  deallocate(M)
  deallocate(u)
  deallocate(v)
  deallocate(p)
  deallocate(scalar_solver)

  call timer(end_time)

  if (par_env%proc_id == par_env%root) then
    print *, "Elapsed time = ", (end_time - start_time)
  end if

  call cleanup_parallel_environment(par_env)

contains

  subroutine initialise_fields(par_env, u, v, p)

    ! Arguments
    class(parallel_environment), intent(in) :: par_env
    class(field), intent(inout) :: u, v, p

    ! Local variables
    integer(accs_int) :: i

    ! Create square mesh
    square_mesh = build_square_mesh(cps, 1.0_accs_real, par_env)

    ! Allocate field arrays
    allocate(u%val(cps,cps))
    allocate(v%val(cps,cps))
    allocate(p%val(cps,cps))

    ! Set initial values for fields
    do i = 1, cps
        u%val(i,:) = real(i, accs_real)/real(cps, accs_real)
        v%val(i,:) = -real(i, accs_real)/real(cps, accs_real)
        p%val(i,:) = 0.0_accs_real
    end do

  end subroutine initialise_fields

end program simple