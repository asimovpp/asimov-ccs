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
       comm=MPI_COMM_NULL, nnz=0)
  type(linear_system) :: poisson_eq
  class(linear_solver), allocatable :: poisson_solver
  
  integer(accs_int), parameter :: cps = 5 ! Cells per side
                                          ! XXX: temporary parameter - this should be read from input

  !!=========================================================
  !! The innards of a rudimentary mesh structure
  integer(accs_int) :: nc
  integer(accs_int) :: istart, iend
  real(accs_real) :: h, Af
  integer(accs_int), dimension(:), allocatable :: nnb
  integer(accs_int), dimension(:, :), allocatable :: nbidx
  !!=========================================================

  real(accs_real) :: err_norm
  
  integer(accs_err) :: ierr
  integer :: comm_rank, comm_size
  
  !! Initialise program
  call accs_init()
  call init_ex3()
  
  !! Create stiffness matrix
  mat_sizes%rglob = nc
  mat_sizes%cglob = mat_sizes%rglob
  mat_sizes%comm = PETSC_COMM_WORLD
  mat_sizes%nnz = 9
  call create_matrix(mat_sizes, M)

  call discretise_poisson(M)
  call update(M) ! Performs the parallel assembly

  !! Create right-hand-side and solution vectors
  vec_sizes%nloc = -1
  vec_sizes%nglob = nc
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

    stop "eval_rhs not converted to FVM"
    
    val_dat%mode = add_mode
    allocate(val_dat%idx(1))
    allocate(val_dat%val(1))
    
    do i = istart, iend
       ii = i - 1
       x = h * modulo(ii, cps)
       y = h * (ii / cps)

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
    
  end subroutine eval_cell_rhs

  subroutine discretise_poisson(M)

    use accs_constants, only : add_mode
    use accs_types, only : matrix_values
    use accs_utils, only : set_values
    
    class(matrix), intent(inout) :: M

    type(matrix_values) :: mat_coeffs
    integer(accs_int) :: i, j

    real(accs_real) :: coeff_f
    
    mat_coeffs%mode = add_mode

    !! Loop over cells
    do i = istart, iend
       allocate(mat_coeffs%rglob(1 + nnb(i)), &
            mat_coeffs%cglob(1 + nnb(i)), &
            mat_coeffs%val(1 + nnb(i)))
       mat_coeffs%rglob(:) = i
       !! Loop over faces
       do j = 1, nnb(i)
          if (nbidx(i, j) >= 0) then
             !! Interior face
             coeff_f = (1.0 / h) * Af
             mat_coeffs%val(1) = mat_coeffs%val(1) - coeff_f
             mat_coeffs%val(j + 1) = coeff_f
             mat_coeffs%cglob(j + 1) = nbidx(i, j)
          end if
       end do
       deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)
       mat_coeffs%cglob = mat_coeffs%rglob
       call set_values(mat_coeffs, M)
    end do

    deallocate(mat_coeffs%rglob)
    deallocate(mat_coeffs%cglob)
    deallocate(mat_coeffs%val)
    
  end subroutine discretise_poisson

  subroutine apply_dirichlet_bcs(M, b, u)

    use accs_constants, only : add_mode
    use accsmat, only : set_eqn
    use accs_types, only : vector_values, matrix_values
    use accs_utils, only : set_values
    
    class(matrix), intent(inout) :: M
    class(vector), intent(inout) :: b, u

    integer(accs_int) :: i, j

    type(vector_values) :: vec_values
    type(matrix_values) :: mat_coeffs

    allocate(mat_coeffs%rglob(1), mat_coeffs%cglob(1), mat_coeffs%val(1))
    allocate(vec_values%idx(1), vec_values%val(1))
    mat_coeffs%mode = add_mode
    vec_values%mode = add_mode
    
    associate(row=>mat_coeffs%rglob, &
         col=>mat_coeffs%cglob, &
         coeff=>mat_coeffs%val, &
         idx=>vec_values%idx, &
         val=>vec_values%val)
      do i = istart, iend
         do j = 1, nnb(i)
            if (nbidx(j, i) < 0) then
               !! Boundary face

               ! Coefficient
               row(1) = i
               col(1) = i
               coeff(1) = -(3 / h) * Af
               call set_values(mat_coeffs, M)
               
               ! RHS vector
               idx(1) = i
               val(1) = h * (i / (cps + 1))
               val(1) = -(val(1) / h) * Af
               call set_values(vec_values, b)
            end if
         end do
      end do
    end associate
    
    deallocate(vec_values%idx)
    deallocate(vec_values%val)

    !! Need to update halo values
    call update(b)
    call update(u)
    
  end subroutine apply_dirichlet_bcs

  subroutine set_exact_sol(ustar)

    use petscvec, only : VecGetOwnershipRange
    
    use accs_constants, only : insert_mode
    use accs_types, only : vector_values
    use accs_utils, only : set_values

    use accs_petsctypes, only : vector_petsc
    
    class(vector), intent(inout) :: ustar

    type(vector_values) :: vec_values
    integer(accs_int) :: i, istart, iend
    integer(accs_err) :: ierr

    stop "set_exact_sol not converted to FVM"
    
    ! TODO: abstract this!
    select type(ustar)
    type is (vector_petsc)
       call VecGetOwnershipRange(ustar%v, istart, iend, ierr)
       istart = istart + 1
    class default
       print *, "Wrong type of vector! Besides you should have abstracted this!"
       stop
    end select
    
    allocate(vec_values%idx(1))
    allocate(vec_values%val(1))
    vec_values%mode = insert_mode
    do i = istart, iend
       vec_values%idx(1) = i - 1
       vec_values%val(1) = h * (vec_values%idx(1) / (cps + 1))
       call set_values(vec_values, ustar)
    end do
    deallocate(vec_values%idx)
    deallocate(vec_values%val)

    call update(ustar)
  end subroutine set_exact_sol

  subroutine init_ex3()

    stop "init_ex3 not converted to FVM"
    
    h = 1.0_accs_real / cps
    nc = cps**2       ! Number of cells
    call MPI_Comm_size(PETSC_COMM_WORLD, comm_size, ierr)
    call MPI_Comm_rank(PETSC_COMM_WORLD, comm_rank, ierr)
    istart = comm_rank * (nc / comm_size)
    if (modulo(nc, comm_size) < comm_rank) then
       istart = istart + modulo(nc, comm_size)
    else
       istart = istart + comm_rank
    end if
    istart = istart + 1 ! Fortran - 1 indexed
  
    iend = istart + nc / comm_size
    if (modulo(nc, comm_size) > comm_rank) then
       iend = iend + 1
    end if
    iend = iend - 1
    
  end subroutine init_ex3
end program ex3
