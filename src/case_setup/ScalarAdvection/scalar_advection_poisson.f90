!> @brief Program file for scalar advection case
!
!

program scalar_advection

  !! ASiMoV-CCS uses
  use kinds, only : accs_real, accs_int
  use types, only : vector_init_data, vector, matrix_init_data, matrix, &
                    linear_system, linear_solver, mesh, set_global_matrix_size, &
                    viewer, field, upwind_field, central_field
  use vec, only : create_vector, axpy, norm!, vec_view
  use mat, only : create_matrix, set_nnz
  use solver, only : create_solver, solve, set_linear_system
  use utils, only : update, begin_update, end_update, finalise, initialise, &
                    set_global_size
  use mesh_utils, only : build_square_mesh
  use petsctypes, only : matrix_petsc
  use parallel_types, only: parallel_environment
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer
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

  real(accs_real) :: err_norm

  double precision :: start_time
  double precision :: end_time

  call initialise_parallel_environment(par_env) 
  call read_command_line_arguments()
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
  call compute_fluxes(M, source, u, v, square_mesh)

  call begin_update(M) ! Start the parallel assembly for M

  call begin_update(source) ! Start the parallel assembly for source
  call end_update(M) ! Complete the parallel assembly for M
  call end_update(source) ! Complete the parallel assembly for source

  !! Create linear solver & set options
  call set_linear_system(scalar_linear_system, source, scalar, M, par_env)
  call create_solver(scalar_linear_system, scalar_solver)
  call solve(scalar_solver)
  
  !call vec_view(scalar)

  !! Check solution
  call set_exact_sol(solution)
  call axpy(-1.0_accs_real, solution, scalar)

  err_norm = norm(scalar, 2) * square_mesh%h
  if (par_env%proc_id == par_env%root) then
     print *, "Norm of error = ", err_norm
  end if
  
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

    !square_mesh = build_square_mesh(cps, 1.0_accs_real, par_env)
    square_mesh = build_square_mesh(cps, 1.0, par_env)

    ! Allocate velocity and scalar field arrays
    allocate(u%val(cps,cps))
    allocate(v%val(cps,cps))

    ! Set IC velocity and scalar fields
    !u%val = 1.0_accs_real
    !v%val = 0.0_accs_real
    do i = 1, cps
      u%val(i,:) = real(i, kind=accs_real)/real(cps, kind=accs_real)
      v%val(:,i) = -real(i, kind=accs_real)/real(cps, kind=accs_real)
    end do
  end subroutine initialise_scalar_advection

  subroutine set_exact_sol(ustar)

    use constants, only : insert_mode
    use types, only : vector_values
    use utils, only : set_values, pack_entries
    
    class(vector), intent(inout) :: ustar

    type(vector_values) :: vec_values
    integer(accs_int) :: i
    
    allocate(vec_values%idx(1))
    allocate(vec_values%val(1))
    vec_values%mode = insert_mode
    do i = 1, square_mesh%nlocal
       associate(idx=>square_mesh%idx_global(i))
         !call pack_entries(vec_values, 1, idx, rhs_val(i))
         call pack_entries(vec_values, 1, idx, 1.0_accs_real)
         call set_values(vec_values, ustar)
       end associate
    end do
    deallocate(vec_values%idx)
    deallocate(vec_values%val)

    call update(ustar)
  end subroutine set_exact_sol

  !> @brief read command line arguments and their values
  subroutine read_command_line_arguments()

    character(len=32) :: arg !> argument string
    integer(accs_int) :: nargs !> number of arguments

    do nargs = 1, command_argument_count()
      call get_command_argument(nargs, arg)
      select case (arg)
        case ('--ccs_m') !> problems size
          call get_command_argument(nargs+1, arg)
          read(arg, '(I5)') cps
        case ('--ccs_help')
          if(par_env%proc_id==par_env%root) then
            print *, "================================"
            print *, "ASiMoV-CCS command line options:"
            print *, "================================"
            print *, "--ccs_help:         This help menu."
            print *, "--ccs_m <value>:    Problem size."
          endif
          stop
        case default
      end select
    end do 

  end subroutine read_command_line_arguments
  
end program scalar_advection
