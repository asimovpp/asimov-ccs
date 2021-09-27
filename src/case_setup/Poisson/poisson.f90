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

  !! External uses
  use MPI
  use petsc, only : PETSC_COMM_WORLD

  !! ASiMoV-CCS uses
  use accs_kinds, only : accs_real, accs_int, accs_err
  use accs_types, only : vector_init_data, vector, matrix_init_data, matrix, linear_system, &
       linear_solver, mesh
  use accsvec, only : create_vector, axpy, norm
  use accsmat, only : create_matrix
  use accs_solver, only : create_solver, solve
  use accs_utils, only : accs_init, accs_finalise, update, finalise

  use accs_mesh, only : build_mesh
  
  use accs_petsctypes, only : matrix_petsc
  
  implicit none

  class(vector), allocatable, target :: u, b
  class(vector), allocatable :: ustar
  type(vector_init_data) :: vec_sizes
  class(matrix), allocatable, target :: M
  type(matrix_init_data) :: mat_sizes = matrix_init_data(rglob=-1, cglob=-1, rloc=-1, cloc=-1, &
       comm=MPI_COMM_NULL, nnz=0)
  type(linear_system) :: poisson_eq
  class(linear_solver), allocatable :: poisson_solver
  type(mesh) :: square_mesh
  
  integer(accs_int), parameter :: cps = 10 ! Cells per side
                                          ! XXX: temporary parameter - this should be read from input

  real(accs_real) :: err_norm

  integer :: comm_rank
  
  !! Initialise program
  call accs_init()
  call init_poisson()

  !! Create stiffness matrix
  mat_sizes%rglob = square_mesh%n
  mat_sizes%cglob = mat_sizes%rglob
  mat_sizes%comm = PETSC_COMM_WORLD
  mat_sizes%nnz = 5
  call create_matrix(mat_sizes, M)

  call discretise_poisson(M)
  call update(M) ! Performs the parallel assembly

  !! Create right-hand-side and solution vectors
  vec_sizes%nloc = -1
  vec_sizes%nglob = square_mesh%n
  vec_sizes%comm = PETSC_COMM_WORLD
  call create_vector(vec_sizes, u)
  call create_vector(vec_sizes, b)
  call update(u) ! Performs the parallel assembly

  !! Evaluate right-hand-side vector
  call eval_rhs(b)
  call update(b) ! Performs the parallel assembly

  !! Modify matrix and right-hand-side vector to apply Dirichlet boundary conditions
  call apply_dirichlet_bcs(M, b)
  call update(b)
  call finalise(M)
  
  !! Create linear solver & set options
  poisson_eq%rhs => b
  poisson_eq%sol => u
  poisson_eq%M => M
  poisson_eq%comm = PETSC_COMM_WORLD
  call create_solver(poisson_eq, poisson_solver)
  call solve(poisson_solver)

  !! Check solution
  call create_vector(vec_sizes, ustar)
  call set_exact_sol(ustar)
  call axpy(-1.0_accs_real, ustar, u)
  err_norm = norm(u, 2) * square_mesh%h
  if (comm_rank == 0) then
     print *, "Norm of error = ", err_norm
  end if
  
  !! Clean up
  deallocate(u)
  deallocate(b)
  deallocate(ustar)
  deallocate(M)
  deallocate(poisson_solver)
  
  call accs_finalise()

contains

  subroutine eval_rhs(b)

    use accs_constants, only : add_mode
    use accs_types, only : vector_values
    use accs_utils, only : set_values
    
    class(vector), intent(inout) :: b

    integer(accs_int) :: i, ii
    real(accs_real) :: x, y, h

    type(vector_values) :: val_dat
    
    val_dat%mode = add_mode
    allocate(val_dat%idx(1))
    allocate(val_dat%val(1))
    
    do i = 1, square_mesh%nlocal
       ii = square_mesh%idx_global(i) - 1              ! This code is translated from C - shift back by 1
       h = square_mesh%h
       x = h * (modulo(ii, cps) + 0.5)
       y = h * ((ii / cps) + 0.5)

       val_dat%idx(1) = ii
       call eval_cell_rhs(x, y, h**2, val_dat%val(1))
       call set_values(val_dat, b)
    end do

    deallocate(val_dat%idx)
    deallocate(val_dat%val)
    
  end subroutine eval_rhs

  !> @brief Apply forcing function
  pure subroutine eval_cell_rhs (x, y, H, r)
    
    real(accs_real), intent(in) :: x, y, H
    real(accs_real), intent(out) :: r
    
    r = 0.0 &
         + 0 * (x + y + H) ! Silence unused dummy argument error
    r = square_mesh%vol * r
    
  end subroutine eval_cell_rhs

  subroutine discretise_poisson(M)

    use accs_constants, only : insert_mode
    use accs_types, only : matrix_values
    use accs_utils, only : set_values
    
    class(matrix), intent(inout) :: M

    type(matrix_values) :: mat_coeffs
    integer(accs_int) :: i, j

    real(accs_real) :: coeff_f
    
    mat_coeffs%mode = insert_mode

    !! Loop over cells
    do i = 1, square_mesh%nlocal
       !> @todo: Doing this in a loop is awful code - malloc maximum coefficients per row once,
       !!        filling from front, and pass the number of coefficients to be set, requires
       !!        modifying the matrix_values type and the implementation of set_values applied to
       !!        matrices.
       associate(idxg=>square_mesh%idx_global(i), &
            nnb=>square_mesh%nnb(i))
         
         allocate(mat_coeffs%rglob(1), &
              mat_coeffs%cglob(1 + nnb))
         allocate(mat_coeffs%val(1 + nnb))

         mat_coeffs%rglob(1) = idxg
         mat_coeffs%cglob(:) = -1 ! -ve indices are ignored
         mat_coeffs%cglob(1) = idxg
         mat_coeffs%val(:) = 0.0_accs_real
       
         !! Loop over faces
         do j = 1, nnb
            associate(nbidxg=>square_mesh%nbidx(j, i))

              if (nbidxg > 0) then
                 !! Interior face
                 coeff_f = (1.0 / square_mesh%h) * square_mesh%Af
                 mat_coeffs%val(1) = mat_coeffs%val(1) - coeff_f
                 mat_coeffs%val(j + 1) = coeff_f
                 mat_coeffs%cglob(j + 1) = nbidxg
              end if

            end associate
         end do
       
         mat_coeffs%rglob(:) = mat_coeffs%rglob(:) - 1
         mat_coeffs%cglob(:) = mat_coeffs%cglob(:) - 1

         call set_values(mat_coeffs, M)

         deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)
         
       end associate

    end do
    
  end subroutine discretise_poisson

  subroutine set_exact_sol(ustar)

    use accs_constants, only : insert_mode
    use accs_types, only : vector_values
    use accs_utils, only : set_values
    
    class(vector), intent(inout) :: ustar

    type(vector_values) :: vec_values
    integer(accs_int) :: i
    
    allocate(vec_values%idx(1))
    allocate(vec_values%val(1))
    vec_values%mode = insert_mode
    do i = 1, square_mesh%nlocal
       associate(idx=>square_mesh%idx_global(i))
         vec_values%idx(1) = idx - 1
         vec_values%val(1) = rhs_val(idx)
         call set_values(vec_values, ustar)
       end associate
    end do
    deallocate(vec_values%idx)
    deallocate(vec_values%val)

    call update(ustar)
  end subroutine set_exact_sol

  subroutine init_poisson()

    integer(accs_err) :: ierr
    
    call MPI_Comm_rank(PETSC_COMM_WORLD, comm_rank, ierr)
    square_mesh = build_mesh(cps, 1.0_accs_real, PETSC_COMM_WORLD)
    
  end subroutine init_poisson

  pure function rhs_val(i, opt_offset) result(r)

    integer(accs_int), intent(in) :: i
    real(accs_real), intent(in), optional :: opt_offset

    integer(accs_int) :: ii
    real(accs_real) :: r
    real(accs_real) :: offset

    if (present(opt_offset)) then
       offset = opt_offset
    else
       offset = 0
    end if

    ii = i - 1

    associate(h=>square_mesh%h)
      r = h * (ii / cps + (0.5 + offset))
    end associate
    
  end function rhs_val
  
end program poisson
