!> @brief Program file PETSc ex3
!>
!> @details Port of PETSc ksp/tutorial/ex3.c to ASiMoV-CCS style code - this is to help
!>          determine how to interface our design with PETSc.

program ex3

  !! External uses
  use MPI
  use petsc, only : PETSC_COMM_WORLD

  !! ASiMoV-CCS uses
  use accs_kinds, only : accs_real, accs_int, accs_err
  use accs_types, only : vector_init_data, vector, matrix_init_data, matrix, linear_system, &
       linear_solver
  use accsvec, only : create_vector, axpy, norm
  use accsmat, only : create_matrix
  use accs_solver, only : create_solver, solve
  use accs_utils, only : accs_init, accs_finalise, update

  use accs_petsctypes, only : matrix_petsc
  
  implicit none

  class(vector), allocatable, target :: u, b
  class(vector), allocatable :: ustar
  type(vector_init_data) :: vec_sizes
  class(matrix), allocatable, target :: M
  type(matrix_init_data) :: mat_sizes = matrix_init_data(rglob=-1, cglob=-1, rloc=-1, cloc=-1, &
       comm=MPI_COMM_NULL)
  type(linear_system) :: poisson_eq
  class(linear_solver), allocatable :: poisson_solver
  
  integer(accs_int), parameter :: eps = 5 ! Elements per side
                                          ! XXX: temporary parameter - this should be read from input
  integer(accs_int) :: nel, nnd
  
  integer(accs_int) :: istart, iend
  real(accs_real) :: h
  integer(accs_int), parameter :: npe = 4 ! Points per element

  real(accs_real) :: err_norm

  real(accs_real), dimension(16) :: K ! Element stiffness matrix
  
  integer(accs_err) :: ierr
  integer :: comm_rank, comm_size
  
  !! Initialise program
  call accs_init()
  h = 1.0_accs_real / eps

  nel = eps**2       ! Number of elements
  nnd = (eps + 1)**2 ! Number of nodes
  call MPI_Comm_size(PETSC_COMM_WORLD, comm_size, ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD, comm_rank, ierr)
  istart = comm_rank * ((nel) / comm_size)
  if (modulo(nel, comm_size) < comm_rank) then
     istart = istart + modulo(nel, comm_size)
  else
     istart = istart + comm_rank
  end if
  istart = istart + 1 ! Fortran - 1 indexed
  
  iend = istart + nel / comm_size
  if (modulo(nel, comm_size) > comm_rank) then
     iend = iend + 1
  end if
  iend = iend - 1
  
  !! Create stiffness matrix
  mat_sizes%rglob = nnd
  mat_sizes%cglob = mat_sizes%rglob
  mat_sizes%comm = PETSC_COMM_WORLD
  call create_matrix(mat_sizes, M)

  call form_element_stiffness(h**2, K)
  call assemble_global_matrix(K, M)
  call update(M) ! Performs the parallel assembly

  !! Create right-hand-side and solution vectors
  vec_sizes%nloc = -1
  vec_sizes%nglob = nnd
  vec_sizes%comm = PETSC_COMM_WORLD
  call create_vector(vec_sizes, u)
  call create_vector(vec_sizes, b)
  call update(u) ! Performs the parallel assembly
  
  !! Evaluate right-hand-side vector
  call eval_rhs(b)
  call update(b) ! Performs the parallel assembly
  
  !! Modify matrix and right-hand-side vector to apply Dirichlet boundary conditions
  call apply_dirichlet_bcs(M, b, u)
  
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
  print *, norm(u, 2) * h, norm(ustar, 2) * h
  call axpy(-1.0_accs_real, ustar, u)
  err_norm = norm(u, 2) * h
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
    real(accs_real) :: x, y

    type(vector_values) :: val_dat

    val_dat%mode = add_mode
    allocate(val_dat%idx(npe))
    allocate(val_dat%val(npe))
    
    do i = istart, iend
       ii = i - 1
       x = h * modulo(ii, eps)
       y = h * (ii / eps)

       call element_indices(i, val_dat%idx)
       call eval_element_rhs(x, y, h**2, val_dat%val)
       call set_values(val_dat, b)
    end do

    deallocate(val_dat%idx)
    deallocate(val_dat%val)
    
  end subroutine eval_rhs

  pure subroutine element_indices (i, idx)

    integer(accs_int), intent(in) :: i
    integer(accs_int), dimension(npe), intent(out) :: idx

    integer(accs_int) :: ii
    
    ii = i - 1 !! Need to convert from Fortran to C indexing as this is based on ex3.c
    idx(1) = (eps + 1) * (ii / eps) + modulo(ii, eps)
    idx(2) = idx(1) + 1
    idx(3) = idx(2) + (eps + 1)
    idx(4) = idx(3) - 1
    
  end subroutine element_indices

  !> @brief Apply forcing function
  pure subroutine eval_element_rhs (x, y, H, r)
    
    real(accs_real), intent(in) :: x, y, H
    real(accs_real), dimension(npe), intent(out) :: r
    
    r(:) = 0.0 &
         + 0 * (x + y + H) ! Silence unused dummy argument error
    
  end subroutine eval_element_rhs

  pure subroutine form_element_stiffness(H, K)

    real(accs_real), intent(in) :: H
    real(accs_real), dimension(16), intent(inout) :: K

    K(1)  = H/6.0;   K(2)  = -.125*H;  K(3)  = H/12.0;  K(4)  = -.125*H;
    K(5)  = -.125*H; K(6)  = H/6.0;    K(7)  = -.125*H; K(8)  = H/12.0;
    K(9)  = H/12.0;  K(10)  = -.125*H; K(11) = H/6.0;   K(12) = -.125*H;
    K(13) = -.125*H; K(14) = H/12.0;   K(15) = -.125*H; K(16) = H/6.0;
    
  end subroutine form_element_stiffness

  subroutine assemble_global_matrix(K, M)

    use accs_constants, only : add_mode
    use accs_types, only : matrix_values
    use accs_utils, only : set_values
    
    real(accs_real), dimension(16), intent(in) :: K
    class(matrix), intent(inout) :: M

    type(matrix_values) :: mat_coeffs
    integer(accs_int) :: i

    allocate(mat_coeffs%rglob(4))
    allocate(mat_coeffs%cglob(4))
    allocate(mat_coeffs%val(16))
    mat_coeffs%val(:) = K(:)
    mat_coeffs%mode = add_mode
    
    do i = istart, iend
       call element_indices(i, mat_coeffs%rglob)
       mat_coeffs%cglob = mat_coeffs%rglob
       call set_values(mat_coeffs, M)
    end do

    deallocate(mat_coeffs%rglob)
    deallocate(mat_coeffs%cglob)
    deallocate(mat_coeffs%val)
  end subroutine assemble_global_matrix

  subroutine apply_dirichlet_bcs(M, b, u)

    use accs_constants, only : insert_mode
    use accsmat, only : set_eqn
    use accs_types, only : vector_values
    use accs_utils, only : set_values
    
    class(matrix), intent(inout) :: M
    class(vector), intent(inout) :: b, u

    integer(accs_int), dimension(4 * eps) :: rows
    integer(accs_int) :: i, idx

    type(vector_values) :: vec_values

    !! Set row indices
    do i = 1, eps + 1
       !! Top of domain
       idx = i
       rows(idx) = i
       !! Bottom of domain
       idx = 3 * eps - 1 + i
       rows(idx) = i + eps * (eps + 1)
    end do
    idx = (eps + 1) + 1
    do i = (eps + 1) + 1, eps * (eps + 1), eps + 1
       rows(idx) = i
       idx = idx + 1
    end do
    idx = 2 * eps + 1
    do i = (2 * eps + 1) + 1, eps * (eps + 1), eps + 1
       rows(idx) = i
       idx = idx + 1
    end do

    !! PETSc is zero-indexed
    rows(:) = rows(:) - 1

    allocate(vec_values%idx(1))
    allocate(vec_values%val(1))
    vec_values%mode = insert_mode
    do i = 1, 4 * eps
       vec_values%idx(1) = rows(i)
       vec_values%val(1) = h * (rows(i) / (eps + 1))
       call set_values(vec_values, b)
       call set_values(vec_values, u)
    end do
    deallocate(vec_values%idx)
    deallocate(vec_values%val)

    call set_eqn(rows, M)

    !! Need to update halo values
    call update(b)
    call update(u)
    
  end subroutine apply_dirichlet_bcs

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
    do i = 1, nnd ! TODO: get the ownership range - this only works for 1 rank!
       vec_values%idx(1) = i - 1
       vec_values%val(1) = h * (vec_values%idx(1) / (eps + 1))
       print *, vec_values%idx(1), vec_values%val(1)
       call set_values(vec_values, ustar)
    end do
    deallocate(vec_values%idx)
    deallocate(vec_values%val)

    call update(ustar)
  end subroutine set_exact_sol
end program ex3
