!> @brief Program file for Poisson case
!
!> @details Based on prototype/ex3 a port of PETSc ksp/tutorial/ex3.c to ASiMoV-CCS style code.
!!          This case demonstrates setting up a linear system and solving it with ASiMoV-CCS, note
!!          the code is independent of PETSc.
!!          The example case solves the equation
!!          \f[
!!            {\nabla^2} p = f
!!          \f]
!!          in the unit square with Dirichlet boundary conditions
!!          \f[
!!            p\left(\boldsymbol{x}\right) = y,\ \boldsymbol{x}\in\partial\Omega
!!          \f]
!

program poisson

  !! ASiMoV-CCS uses
  use kinds, only : accs_real, accs_int
  use types, only : vector_init_data, vector, matrix_init_data, matrix, &
                    linear_system, linear_solver, mesh, set_global_matrix_size
  use vec, only : create_vector, axpy, norm
  use mat, only : create_matrix, set_nnz
  use solver, only : create_solver, solve, set_linear_system
  use utils, only : update, begin_update, end_update, finalise, initialise, &
                    set_global_size
  use mesh_utils, only : build_square_mesh
  use petsctypes, only : matrix_petsc
  use parallel_types, only: parallel_environment
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer
  
  implicit none

  class(parallel_environment), allocatable, target :: par_env
  class(vector), allocatable, target :: u, b
  class(vector), allocatable :: ustar
  class(matrix), allocatable, target :: M
  class(linear_solver), allocatable :: poisson_solver

  type(vector_init_data) :: vec_sizes
  type(matrix_init_data) :: mat_sizes
  type(linear_system) :: poisson_eq
  type(mesh) :: square_mesh

  integer(accs_int) :: cps = 10 ! Default value for cells per side

  real(accs_real) :: err_norm

  double precision :: start_time
  double precision :: end_time

  call initialise_parallel_environment(par_env) 
  call read_command_line_arguments()
  call timer(start_time)

  call initialise_poisson(par_env)

  !! Initialise with default values
  call initialise(mat_sizes)
  call initialise(vec_sizes)
  call initialise(poisson_eq)

  !! Create stiffness matrix
  call set_global_size(mat_sizes, square_mesh%n, square_mesh%n, par_env)
  call set_nnz(mat_sizes, 5)
  call create_matrix(mat_sizes, M)

  call discretise_poisson(M)

  call begin_update(M) ! Start the parallel assembly for M

  !! Create right-hand-side and solution vectors
  call set_global_size(vec_sizes, square_mesh%n, par_env)
  call create_vector(vec_sizes, b)
  call create_vector(vec_sizes, ustar)
  call create_vector(vec_sizes, u)

  call begin_update(u) ! Start the parallel assembly for u

  !! Evaluate right-hand-side vector
  call eval_rhs(b)

  call begin_update(b) ! Start the parallel assembly for b
  call end_update(M) ! Complete the parallel assembly for M
  call end_update(b) ! Complete the parallel assembly for b

  !! Modify matrix and right-hand-side vector to apply Dirichlet boundary conditions
  call apply_dirichlet_bcs(M, b)
  call begin_update(b) ! Start the parallel assembly for b
  call finalise(M)

  call end_update(u) ! Complete the parallel assembly for u
  call end_update(b) ! Complete the parallel assembly for b

  !! Create linear solver & set options
  call set_linear_system(poisson_eq, b, u, M, par_env)
  call create_solver(poisson_eq, poisson_solver)
  call solve(poisson_solver)

  !! Check solution
  call set_exact_sol(ustar)
  call axpy(-1.0_accs_real, ustar, u)

  err_norm = norm(u, 2) * square_mesh%h
  if (par_env%proc_id == par_env%root) then
     print *, "Norm of error = ", err_norm
  end if
  
  !! Clean up
  deallocate(u)
  deallocate(b)
  deallocate(ustar)
  deallocate(M)
  deallocate(poisson_solver)

  call timer(end_time)

  if (par_env%proc_id == par_env%root) then
    print *, "Elapsed time = ", (end_time - start_time)
  end if

  call cleanup_parallel_environment(par_env)

contains

  subroutine eval_rhs(b)

    use constants, only : add_mode
    use types, only : vector_values
    use utils, only : set_values, pack_entries
    
    class(vector), intent(inout) :: b

    integer(accs_int) :: i
    integer(accs_int) :: nloc
    real(accs_real) :: h
    real(accs_real) :: r

    type(vector_values) :: val_dat
    
    val_dat%mode = add_mode
    allocate(val_dat%idx(1))
    allocate(val_dat%val(1))

    nloc = square_mesh%nlocal

    h = square_mesh%h

    ! this is currently setting 1 vector value at a time
    ! consider changing to doing all the updates in one go
    ! to do only 1 call to eval_cell_rhs and set_values
    associate(idx => square_mesh%idx_global)
      do i = 1, nloc
        associate(x => square_mesh%xc(1, i), &
             y => square_mesh%xc(2, i), &
             V => square_mesh%vol(i))
          call eval_cell_rhs(x, y, h**2, r)
          r = V * r
          call pack_entries(val_dat, 1, idx(i), r)
          call set_values(val_dat, b)
        end associate
      end do
    end associate
    
    deallocate(val_dat%idx)
    deallocate(val_dat%val)
    
  end subroutine eval_rhs

  !> @brief Apply forcing function
  pure subroutine eval_cell_rhs (x, y, H, r)
    
    real(accs_real), intent(in) :: x, y, H
    real(accs_real), intent(out) :: r
    
    r = 0.0 &
         + 0 * (x + y + H) ! Silence unused dummy argument error
    
  end subroutine eval_cell_rhs

  subroutine discretise_poisson(M)

    use constants, only : insert_mode
    use types, only : matrix_values
    use utils, only : set_values, pack_entries
    
    class(matrix), intent(inout) :: M

    type(matrix_values) :: mat_coeffs
    integer(accs_int) :: i, j

    integer(accs_int) :: row, col
    real(accs_real) :: coeff_f, coeff_p, coeff_nb
    
    mat_coeffs%mode = insert_mode

    !! Loop over cells
    do i = 1, square_mesh%nlocal
      !> @todo: Doing this in a loop is awful code - malloc maximum coefficients per row once,
      !!        filling from front, and pass the number of coefficients to be set, requires
      !!        modifying the matrix_values type and the implementation of set_values applied to
      !!        matrices.
      associate(idxg=>square_mesh%idx_global(i), &
                nnb=>square_mesh%nnb(i))
        
        allocate(mat_coeffs%rglob(1))
        allocate(mat_coeffs%cglob(1 + nnb))
        allocate(mat_coeffs%val(1 + nnb))

        row = idxg
        coeff_p = 0.0_accs_real
      
        !! Loop over faces
        do j = 1, nnb
          coeff_f = (1.0 / square_mesh%h) * square_mesh%Af(j, i)
          associate(nbidxg=>square_mesh%nbidx(j, i))

            if (nbidxg > 0) then
              !! Interior face
              coeff_p = coeff_p - coeff_f
              coeff_nb = coeff_f
              col = nbidxg
            else
              col = -1
              coeff_nb = 0.0_accs_real
            end if
            call pack_entries(mat_coeffs, 1, j + 1, row, col, coeff_nb)

          end associate
        end do

        !! Add the diagonal entry
        col = row
        call pack_entries(mat_coeffs, 1, 1, row, col, coeff_p)

        !! Set the values
        call set_values(mat_coeffs, M)

        deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)
        
      end associate

    end do
    
  end subroutine discretise_poisson

  subroutine apply_dirichlet_bcs(M, b)

    use constants, only : add_mode
    use mat, only : set_eqn
    use types, only : vector_values, matrix_values, matrix, vector, mesh
    use utils, only : set_values, pack_entries
    use kinds, only: accs_int, accs_real
  
    implicit none
  
    class(matrix), intent(inout) :: M
    class(vector), intent(inout) :: b
  
    integer(accs_int) :: i, j
    real(accs_real) :: boundary_coeff, boundary_val

    integer(accs_int) :: idx, row, col
    real(accs_real) :: r, coeff
    
    type(vector_values) :: vec_values
    type(matrix_values) :: mat_coeffs
  
    allocate(mat_coeffs%rglob(1))
    allocate(mat_coeffs%cglob(1))
    allocate(mat_coeffs%val(1))
    allocate(vec_values%idx(1))
    allocate(vec_values%val(1))

    mat_coeffs%mode = add_mode
    vec_values%mode = add_mode

    associate(idx_global=>square_mesh%idx_global)
  
      do i = 1, square_mesh%nlocal
        if (minval(square_mesh%nbidx(:, i)) < 0) then
          coeff = 0.0_accs_real 
          r = 0.0_accs_real
          
          row = idx_global(i)
          col = idx_global(i)
          idx = idx_global(i)

          do j = 1, square_mesh%nnb(i)

            associate(nbidx=>square_mesh%nbidx(j, i))

              if (nbidx < 0) then
                boundary_coeff = (2.0 / square_mesh%h) * square_mesh%Af(j, i)
                boundary_val = rhs_val(i, j)

                ! Coefficient
                coeff = coeff - boundary_coeff

                ! RHS vector
                r = r - boundary_val * boundary_coeff
              end if

            end associate
          end do

          call pack_entries(mat_coeffs, 1, 1, row, col, coeff)
          call pack_entries(vec_values, 1, idx, r)

          call set_values(mat_coeffs, M)
          call set_values(vec_values, b)
          
        end if
      end do
  
    end associate
  
    deallocate(vec_values%idx)
    deallocate(vec_values%val)
  
  end subroutine apply_dirichlet_bcs

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
         call pack_entries(vec_values, 1, idx, rhs_val(i))
         call set_values(vec_values, ustar)
       end associate
    end do
    deallocate(vec_values%idx)
    deallocate(vec_values%val)

    call update(ustar)
  end subroutine set_exact_sol

  subroutine initialise_poisson(par_env)

    class(parallel_environment) :: par_env

    square_mesh = build_square_mesh(cps, 1.0_accs_real, par_env)
    
  end subroutine initialise_poisson

  pure function rhs_val(i, f) result(r)

    integer(accs_int), intent(in) :: i !> Cell index
    integer(accs_int), intent(in), optional :: f !> Face index (local wrt cell)

    real(accs_real) :: r

    if (present(f)) then
      !! Face-centred value
      associate(y => square_mesh%xf(2, f, i))
        r = y
      end associate
    else
      !! Cell-centred value
      associate(y => square_mesh%xc(2, i))
        r = y
      end associate
    end if
    
  end function rhs_val

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
  
end program poisson
