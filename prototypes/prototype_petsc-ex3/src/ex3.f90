!> @brief Program file PETSc ex3
!
!> @details Port of PETSc ksp/tutorial/ex3.c to ASiMoV-CCS style code - this is to help
!!          determine how to interface our design with PETSc.
!!          Solves the equation
!!          \f[
!!            {\nabla^2} p = f
!!          \f]
!!          in the unit square with Dirichlet boundary conditions
!!          \f[
!!            p\left(\boldsymbol{x}\right) = y,\ \boldsymbol{x}\in\partial\Omega
!!          \f]

program ex3

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
  type(mesh) :: square_mesh
  
  integer(accs_int), parameter :: cps = 5 ! Cells per side
                                          ! XXX: temporary parameter - this should be read from input

  real(accs_real) :: err_norm
  
  integer(accs_err) :: ierr
  integer :: comm_rank, comm_size
  
  !! Initialise program
  call accs_init()
  call init_ex3()

  !! Create stiffness matrix
  mat_sizes%rglob = square_mesh%n
  mat_sizes%cglob = mat_sizes%rglob
  mat_sizes%comm = PETSC_COMM_WORLD
  mat_sizes%nnz = 5
  call create_matrix(mat_sizes, M)

  call discretise_poisson(M)
  print *, "Done"
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
  call apply_dirichlet_bcs(M, b, u)
  call update(M)
  
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
       allocate(mat_coeffs%rglob(1), &
            mat_coeffs%cglob(1 + square_mesh%nnb(i)))
       allocate(mat_coeffs%val(1 + square_mesh%nnb(i)))
       
       mat_coeffs%rglob(1) = square_mesh%idx_global(i)
       mat_coeffs%cglob(:) = -1 ! -ve indices are ignored
       mat_coeffs%cglob(1) = square_mesh%idx_global(i)
       mat_coeffs%val(:) = 0.0_accs_real
       !! Loop over faces
       do j = 1, square_mesh%nnb(i)
          if (square_mesh%nbidx(j, i) >= 0) then
             !! Interior face
             coeff_f = (1.0 / square_mesh%h) * square_mesh%Af
             mat_coeffs%val(1) = mat_coeffs%val(1) - coeff_f
             mat_coeffs%val(j + 1) = coeff_f
             mat_coeffs%cglob(j + 1) = square_mesh%nbidx(j, i)
          end if
       end do
       mat_coeffs%rglob(:) = mat_coeffs%rglob(:) - 1
       mat_coeffs%cglob(:) = mat_coeffs%cglob(:) - 1
       call set_values(mat_coeffs, M)

       deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)

    end do
    print *, "Finished disc Poiss"
    
  end subroutine discretise_poisson

  subroutine apply_dirichlet_bcs(M, b, u)

    use accs_constants, only : add_mode
    use accsmat, only : set_eqn
    use accs_types, only : vector_values, matrix_values
    use accs_utils, only : set_values
    
    class(matrix), intent(inout) :: M
    class(vector), intent(inout) :: b, u

    integer(accs_int) :: i, j
    real(accs_real) :: h

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
      do i = 1, square_mesh%nlocal
         do j = 1, square_mesh%nnb(i)
            if (square_mesh%nbidx(j, i) < 0) then
               !! Boundary face
               h = square_mesh%h

               ! Coefficient
               row(1) = i - 1
               col(1) = i - 1
               coeff(1) = -(3 / h) * square_mesh%Af
               call set_values(mat_coeffs, M)
               
               ! RHS vector
               idx(1) = i - 1
               val(1) = h * (i / (cps + 1) + 0.5) ! XXX: This is sort of y (y=i/cps+0.5), but matches ex3.c...
               val(1) = -(val(1) / h) * square_mesh%Af
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
    
    !> @todo: abstract this!
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
       vec_values%val(1) = square_mesh%h * (vec_values%idx(1) / (cps + 1) + 0.5)
       call set_values(vec_values, ustar)
    end do
    deallocate(vec_values%idx)
    deallocate(vec_values%val)

    call update(ustar)
  end subroutine set_exact_sol

  subroutine init_ex3()

    square_mesh = build_mesh(cps**2)
    
  end subroutine init_ex3

  function build_mesh(n) result(square_mesh)

    integer(accs_int), intent(in) :: n

    integer(accs_int) :: istart, iend
    integer(accs_int) :: i, ii
    type(mesh) :: square_mesh

    square_mesh%n = n             ! Number of cells
    square_mesh%h = 1.0_accs_real / sqrt(real(n, accs_real))
    square_mesh%Af = square_mesh%h     ! Face area
    square_mesh%vol = square_mesh%h**2 ! Cell volume

    !! Setup ownership range
    call MPI_Comm_size(PETSC_COMM_WORLD, comm_size, ierr)
    call MPI_Comm_rank(PETSC_COMM_WORLD, comm_rank, ierr)
    istart = comm_rank * (n / comm_size)
    if (modulo(square_mesh%n, comm_size) < comm_rank) then
       istart = istart + modulo(n, comm_size)
    else
       istart = istart + comm_rank
    end if
    istart = istart + 1 ! Fortran - 1 indexed
  
    iend = istart + n / comm_size
    if (modulo(square_mesh%n, comm_size) > comm_rank) then
       iend = iend + 1
    end if
    iend = iend - 1

    square_mesh%nlocal = (iend - istart) + 1
    allocate(square_mesh%idx_global(square_mesh%nlocal), square_mesh%nnb(square_mesh%nlocal))
    square_mesh%nnb(:) = 4                    ! All cells have 4 neighbours (possibly ghost/boundary cells)
    allocate(square_mesh%nbidx(4, square_mesh%nlocal))
    
    !! Get neighbour indices
    !! XXX: These are global indices and thus may be off-process
    do i = istart, iend
       square_mesh%idx_global(i) = i
       ii = i - 1
       
       !! Left neighbour
       if (modulo(ii, cps) == 0) then
          square_mesh%nbidx(1, i) = -1
       else
          square_mesh%nbidx(1, i) = i - 1
       end if

       !! Right neighbour
       if (modulo(ii, cps) == (cps - 1)) then
          square_mesh%nbidx(2, i) = -2
       else
          square_mesh%nbidx(2, i) = i + 1
       end if

       !! Down neighbour
       if ((ii / cps) == 0) then
          square_mesh%nbidx(3, i) = -3
       else
          square_mesh%nbidx(3, i) = i - cps
       end if

       !! Up neighbour
       if ((ii / cps) == (cps - 1)) then
          square_mesh%nbidx(4, i) = -4
       else
          square_mesh%nbidx(4, i) = i + cps
       end if
    end do
    
  end function build_mesh
  
end program ex3
