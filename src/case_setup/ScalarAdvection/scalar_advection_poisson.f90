!> @brief Program file for scalar advection case
!
!

program scalar_advection

  !! ASiMoV-CCS uses
  use kinds, only : accs_real, accs_int
  use types, only : vector_init_data, vector, matrix_init_data, matrix, &
                    linear_system, linear_solver, mesh, set_global_matrix_size, viewer
  use vec, only : create_vector, axpy, norm, vec_view
  use mat, only : create_matrix, set_nnz
  use solver, only : create_solver, solve, set_linear_system
  use utils, only : update, begin_update, end_update, finalise, initialise, &
                    set_global_size
  use mesh_utils, only : build_square_mesh
  use petsctypes, only : matrix_petsc
  use parallel_types, only: parallel_environment
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer

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

  real(accs_real), dimension(:,:), allocatable :: u, v          ! Prescribed x, y velocity fields

  integer(accs_int) :: cps = 50 ! Default value for cells per side

  real(accs_real) :: err_norm

  double precision :: start_time
  double precision :: end_time

  call initialise_parallel_environment(par_env) 
  call read_command_line_arguments()
  call timer(start_time)

  ! Init ICs (velocities, BC scalar, mesh, etc)
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
  call compute_fluxes(M, source, u, v, 0.0_accs_real, 1.0_accs_real)

  call begin_update(M) ! Start the parallel assembly for M

  call begin_update(scalar) ! Start the parallel assembly for scalar

  call begin_update(source) ! Start the parallel assembly for source
  call end_update(M) ! Complete the parallel assembly for M
  call end_update(source) ! Complete the parallel assembly for source
  call end_update(scalar) ! Complete the parallel assembly for scalar

  !! Create linear solver & set options
  call set_linear_system(scalar_linear_system, source, scalar, M, par_env)
  call create_solver(scalar_linear_system, scalar_solver)
  call solve(scalar_solver)
  
  call vec_view(scalar)

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

  ! ALEXEI: working implementation, but trying to tidy it up
  !subroutine compute_fluxes(M, b, u, v, n_value, w_value)
  !  use constants, only : insert_mode, add_mode
  !  use types, only : matrix_values, vector_values
  !  use utils, only : set_values, pack_entries

  !  class(matrix), intent(inout) :: M   
  !  class(vector), intent(inout) :: b   
  !  real(accs_real), dimension(:,:), intent(in) :: u, v
  !  real(accs_real), intent(in) :: n_value, w_value
  !  
  !  integer(accs_int) :: j           ! Counters
  !  real(accs_real) :: adv_coeff, diff_coeff
  !  real(accs_real) :: adv_coeff_total, diff_coeff_total
  !  integer(accs_int) :: ngb_idx, self_idx
  !  integer(accs_int) :: cps, row, col
  !  real(accs_real) :: BC_value
  !  real(accs_real) :: face_area

  !  type(matrix_values) :: mat_coeffs
  !  type(vector_values) :: b_coeffs
  !  
  !  mat_coeffs%mode = add_mode
  !  b_coeffs%mode = add_mode

  !  allocate(mat_coeffs%rglob(1))
  !  allocate(mat_coeffs%cglob(1))
  !  allocate(mat_coeffs%val(1))
  !  allocate(b_coeffs%idx(1))
  !  allocate(b_coeffs%val(1))

  !  cps = int(sqrt(real(square_mesh%n)))

  !  ! Loop over cells computing advection and diffusion fluxes
  !  do self_idx = 1, square_mesh%n
  !    ! Calculate contribution from neighbours
  !    diff_coeff_total = 0.0_accs_real
  !    adv_coeff_total = 0.0_accs_real
  !    do j = 1, square_mesh%nnb(self_idx)
  !      ngb_idx = square_mesh%nbidx(j, self_idx)
  !      face_area = square_mesh%Af(j, self_idx)
  !      call calc_diffusion_coeff(diff_coeff)
  !      ! Only setting matrix elements that correspond to entries within the domain is the same as zero-gradient BCs when using a uniform square mesh.
  !      if (ngb_idx > 0) then
  !        call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, "CDS", cps, u, v, 0)
  !        call pack_entries(mat_coeffs, 1, 1, self_idx, ngb_idx, adv_coeff + diff_coeff)
  !        call set_values(mat_coeffs, M)
  !        diff_coeff_total = diff_coeff_total + diff_coeff
  !        adv_coeff_total = adv_coeff_total + adv_coeff
  !      else if (ngb_idx == -1) then
  !        ! Set Dirichlet BCs at left boundary
  !        call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, "CDS", cps, u, v, ngb_idx)
  !        call calc_cell_coords(self_idx, cps, row, col)
  !        BC_value = -(1.0_accs_real - real(row, kind=accs_real)/real(cps, kind=accs_real)) * w_value
  !        call pack_entries(b_coeffs, 1, self_idx, (adv_coeff + diff_coeff)*BC_value)
  !        call set_values(b_coeffs, b)
  !        diff_coeff_total = diff_coeff_total + diff_coeff
  !        adv_coeff_total = adv_coeff_total + adv_coeff
  !      else if (ngb_idx == -4) then
  !        ! Set Dirichlet BCs at top boundary
  !        call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, "CDS", cps, u, v, ngb_idx)
  !        call pack_entries(b_coeffs, 1, self_idx, (adv_coeff + diff_coeff)*n_value)
  !        call set_values(b_coeffs, b)
  !        diff_coeff_total = diff_coeff_total + diff_coeff
  !        adv_coeff_total = adv_coeff_total + adv_coeff
  !      end if
  !    end do
  !    call pack_entries(mat_coeffs, 1, 1, self_idx, self_idx, -(adv_coeff_total + diff_coeff_total))
  !    call set_values(mat_coeffs, M)
  !  end do

  !  deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)
  !  deallocate(b_coeffs%idx, b_coeffs%val)
  !end subroutine compute_fluxes

  subroutine compute_fluxes(M, b, u, v, n_value, w_value)
    use constants, only : insert_mode, add_mode
    use types, only : matrix_values, vector_values
    use utils, only : set_values, pack_entries

    class(matrix), intent(inout) :: M   
    class(vector), intent(inout) :: b   
    real(accs_real), dimension(:,:), intent(in) :: u, v
    real(accs_real), intent(in) :: n_value, w_value
    
    integer(accs_int) :: j           ! Counters
    real(accs_real) :: adv_coeff, diff_coeff
    real(accs_real) :: adv_coeff_total, diff_coeff_total
    integer(accs_int) :: ngb_idx, self_idx
    integer(accs_int) :: cps, row, col
    real(accs_real) :: BC_value
    real(accs_real) :: face_area
    integer(accs_int) :: assigned
    integer(accs_int), dimension(:,:), allocatable :: BC_ids
    integer(accs_int) :: b_nnz, mat_nnz

    type(matrix_values) :: mat_coeffs, mat_coeff
    type(matrix_values) :: diagonal_coeffs
    type(vector_values) :: b_coeffs

    integer(accs_int) :: n_bc_cells, n_int_cells
    integer(accs_int) :: bc_counter, mat_counter
    
    cps = int(sqrt(real(square_mesh%n)))

    !allocate(BC_ids(4*cps, 2))
    
    mat_coeffs%mode = add_mode
    mat_coeff%mode = add_mode
    diagonal_coeffs%mode = add_mode
    b_coeffs%mode = add_mode

    mat_nnz = calc_matrix_nnz(cps)
    b_nnz = calc_rhs_nnz(cps)

    !allocate(mat_coeffs%rglob(1))
    !allocate(mat_coeffs%cglob(1 + mat_nnz))
    !allocate(mat_coeffs%val(1 + mat_nnz))
    !mat_coeffs%rglob(:) = 1
    !mat_coeffs%cglob(:) = 1
    !mat_coeffs%val(:) = 0.0_accs_real
    !allocate(diagonal_coeffs%rglob(1))
    !allocate(diagonal_coeffs%cglob(square_mesh%n))
    !allocate(diagonal_coeffs%val(square_mesh%n))
    !diagonal_coeffs%rglob(:) = 1
    !diagonal_coeffs%cglob(:) = 1
    !diagonal_coeffs%val(:) = 0.0_accs_real
    !allocate(b_coeffs%idx(b_nnz))
    !allocate(b_coeffs%val(b_nnz))
    !b_coeffs%idx(:) = 1
    !b_coeffs%val(:) = 0.0_accs_real
    !allocate(mat_coeffs%rglob(1))
    !allocate(mat_coeffs%cglob(1))
    !allocate(mat_coeffs%val(1))
    !allocate(mat_coeff%rglob(1))
    !allocate(mat_coeff%cglob(1))
    !allocate(mat_coeff%val(1))
    !allocate(b_coeffs%idx(1))
    !allocate(b_coeffs%val(1))


    ! Loop over cells computing advection and diffusion fluxes
    !call compute_interior_coeffs(mat_coeffs, diagonal_coeffs, BC_ids, "CDS", u, v)
    n_int_cells = calc_matrix_nnz(cps)
    call compute_interior_coeffs(M, "CDS", u, v, n_int_cells)

    ! Loop over boundaries
    !call compute_boundary_coeffs(mat_coeffs, diagonal_coeffs, b_coeffs, BC_ids, "CDS", u, v, cps)
    n_bc_cells = calc_rhs_nnz(cps)
    call compute_boundary_coeffs(M, b, "CDS", u, v, n_bc_cells)

    !! ALEXEI: Debugging
    !!do j = 1,mat_nnz
    !!  print *, 'j mat val', j, mat_coeffs%cglob(j), mat_coeffs%rglob(j)
    !!end do
    !!stop

    !! Assign values to matrix
    !print *, 'here 1'
    !call set_values(mat_coeffs, M)
    !print *, 'here 2'
    !call set_values(diagonal_coeffs, M)
    !print *, 'here 3'
    !call set_values(b_coeffs, b)
    !print *, 'here 4'


    !bc_counter = 1
    !do self_idx = 1, square_mesh%n
    !  mat_counter = 1
    !  ! Calculate contribution from neighbours
    !  diff_coeff_total = 0.0_accs_real
    !  adv_coeff_total = 0.0_accs_real
    !  do j = 1, square_mesh%nnb(self_idx)
    !    ngb_idx = square_mesh%nbidx(j, self_idx)
    !    face_area = square_mesh%Af(j, self_idx)
    !    call calc_diffusion_coeff(diff_coeff)
    !    ! Only setting matrix elements that correspond to entries within the domain is the same as zero-gradient BCs when using a uniform square mesh.
    !    if (ngb_idx > 0) then
    !      call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, "CDS", cps, u, v, 0)
    !      call pack_entries(mat_coeffs, 1, mat_counter, self_idx, ngb_idx, adv_coeff + diff_coeff)
    !      !call set_values(mat_coeff, M)
    !      !print *, 'mat_counter ', mat_counter, ' row col ', self_idx, ngb_idx
    !      mat_counter = mat_counter + 1
    !      diff_coeff_total = diff_coeff_total + diff_coeff
    !      adv_coeff_total = adv_coeff_total + adv_coeff
    !    else if (ngb_idx == -1) then
    !      ! Set Dirichlet BCs at left boundary
    !      call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, "CDS", cps, u, v, ngb_idx)
    !      call calc_cell_coords(self_idx, cps, row, col)
    !      BC_value = -(1.0_accs_real - real(row, kind=accs_real)/real(cps, kind=accs_real)) * w_value
    !      call pack_entries(b_coeffs, bc_counter, self_idx, (adv_coeff + diff_coeff)*BC_value)
    !      !print *, 'bc_counter self_idx ', bc_counter, self_idx
    !      bc_counter = bc_counter + 1
    !      diff_coeff_total = diff_coeff_total + diff_coeff
    !      adv_coeff_total = adv_coeff_total + adv_coeff
    !    else if (ngb_idx == -4) then
    !      ! Set Dirichlet BCs at top boundary
    !      call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, "CDS", cps, u, v, ngb_idx)
    !      call pack_entries(b_coeffs, bc_counter, self_idx, (adv_coeff + diff_coeff)*n_value)
    !      !print *, 'bc_counter self_idx ', bc_counter, self_idx
    !      bc_counter = bc_counter + 1
    !      diff_coeff_total = diff_coeff_total + diff_coeff
    !      adv_coeff_total = adv_coeff_total + adv_coeff
    !    end if
    !  end do
    !  call pack_entries(mat_coeff, 1, 1, self_idx, self_idx, -(adv_coeff_total + diff_coeff_total))
    !  call set_values(mat_coeff, M)
    !  !call pack_entries(mat_coeffs, mat_counter, mat_counter, self_idx, self_idx, -(adv_coeff_total + diff_coeff_total))
    !  !print *, 'mat_counter ', mat_counter, ' row col ', self_idx, self_idx
    !  !mat_counter = mat_counter + 1
    !  call set_values(mat_coeffs, M)
    !end do
    !call set_values(b_coeffs, b)
    !print *, 'mat_counter ', mat_counter

    !deallocate(BC_ids)
    !deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)
    !deallocate(mat_coeff%rglob, mat_coeff%cglob, mat_coeff%val)
    !!deallocate(diagonal_coeffs%rglob, diagonal_coeffs%cglob, diagonal_coeffs%val)
    !deallocate(b_coeffs%idx, b_coeffs%val)
  end subroutine compute_fluxes

  ! Note: this assumes a 2d grid
  pure function calc_matrix_nnz(cps) result(nnz)
    integer(accs_int), intent(in) :: cps
    integer(accs_int) :: nnz

    nnz = 0

    ! First count the interior cells, each cell has 4 neighbours so contributes 5 entries
    nnz = nnz + 5*(cps-2)*(cps-2)

    ! Add edge cells (excludes corners), each cell has 3 neighbours so contributes 4 entries, 4 edges
    nnz = nnz + 16*(cps-2)

    ! Finally add corner cells, each cell has 2 neighbours so contributes 3 entries, 4 corners
    nnz = nnz + 12
    nnz = 5
  end function calc_matrix_nnz

  ! Note: this assumes a 2d grid
  pure function calc_rhs_nnz(cps) result(nnz)
    implicit none
    integer(accs_int), intent(in) :: cps
    integer(accs_int) :: nnz

    nnz = 2*cps
  end function calc_rhs_nnz

  subroutine compute_interior_coeffs(mat, discretisation, u, v, n_int_cells)
    use constants, only : insert_mode, add_mode
    use types, only: matrix_values
    use utils, only: pack_entries, set_values

    class(matrix), intent(inout) :: mat
    character(len=3), intent(in) :: discretisation
    real(accs_real), dimension(:,:), intent(in) :: u, v
    integer(accs_int), intent(in) :: n_int_cells

    type(matrix_values) :: mat_coeffs
    integer(accs_int) :: self_idx, ngb_idx
    integer(accs_int) :: j
    real(accs_real) :: face_area
    real(accs_real) :: diff_coeff, diff_coeff_total
    real(accs_real) :: adv_coeff, adv_coeff_total
    integer(accs_int) :: mat_counter

    mat_coeffs%mode = add_mode

    allocate(mat_coeffs%rglob(1))
    allocate(mat_coeffs%cglob(n_int_cells))
    allocate(mat_coeffs%val(n_int_cells))

    do self_idx = 1, square_mesh%n
      ! Calculate contribution from neighbours
      mat_counter = 1
      adv_coeff_total = 0.0_accs_real
      diff_coeff_total = 0.0_accs_real
      do j = 1, square_mesh%nnb(self_idx)
        ngb_idx = square_mesh%nbidx(j, self_idx)
        face_area = square_mesh%Af(j, self_idx)
        call calc_diffusion_coeff(diff_coeff)
        if (ngb_idx > 0) then
          call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, discretisation, cps, u, v, 0)
          call pack_entries(mat_coeffs, 1, mat_counter, self_idx, ngb_idx, adv_coeff + diff_coeff)
          !call set_values(mat_coeffs, mat)
          mat_counter = mat_counter + 1
          adv_coeff_total = adv_coeff_total + adv_coeff
          diff_coeff_total = diff_coeff_total + diff_coeff
        else
          call pack_entries(mat_coeffs, 1, mat_counter, self_idx, -1, 0.0_accs_real)
          mat_counter = mat_counter + 1
        end if
      end do
      call pack_entries(mat_coeffs, 1, mat_counter, self_idx, self_idx, -(adv_coeff_total + diff_coeff_total))
      !call set_values(mat_coeffs, mat)
      mat_counter = mat_counter + 1
      call set_values(mat_coeffs, mat)
    end do
    print *, 'mat_counter ', mat_counter

    deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)
  end subroutine compute_interior_coeffs

  subroutine compute_boundary_coeffs(M, b, discretisation, u, v, n_bc_cells)
    use constants, only : insert_mode, add_mode
    use types, only: matrix_values, vector_values
    use utils, only: pack_entries, set_values

    class(matrix), intent(inout) :: M
    class(vector), intent(inout) :: b
    character(len=3), intent(in) :: discretisation
    real(accs_real), dimension(:,:), intent(in) :: u, v
    integer(accs_int), intent(in) :: n_bc_cells

    type(matrix_values) :: mat_coeffs
    type(vector_values) :: b_coeffs

    integer(accs_int) :: self_idx, ngb_idx
    integer(accs_int) :: j
    integer(accs_int) :: bc_counter
    integer(accs_int) :: row, col
    real(accs_real) :: face_area
    real(accs_real) :: diff_coeff
    real(accs_real) :: adv_coeff
    real(accs_real) :: BC_value
    real(accs_real) :: n_value, w_value

    mat_coeffs%mode = add_mode
    b_coeffs%mode = add_mode

    allocate(mat_coeffs%rglob(n_bc_cells))
    allocate(mat_coeffs%cglob(1))
    allocate(mat_coeffs%val(n_bc_cells))
    allocate(b_coeffs%idx(n_bc_cells))
    allocate(b_coeffs%val(n_bc_cells))

    n_value = 0.0_accs_real
    w_value = 1.0_accs_real
    bc_counter = 1
    call calc_diffusion_coeff(diff_coeff)
    do self_idx = 1, square_mesh%n
      ! Calculate contribution from neighbours
      do j = 1, square_mesh%nnb(self_idx)
        ngb_idx = square_mesh%nbidx(j, self_idx)
        face_area = square_mesh%Af(j, self_idx)
        if (ngb_idx == -1) then
          ! Set Dirichlet BCs at left boundary
          call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, "CDS", cps, u, v, ngb_idx)
          call calc_cell_coords(self_idx, cps, row, col)
          BC_value = -(1.0_accs_real - real(row, kind=accs_real)/real(cps, kind=accs_real)) * w_value
          call pack_entries(b_coeffs, bc_counter, self_idx, (adv_coeff + diff_coeff)*BC_value)
          call pack_entries(mat_coeffs, bc_counter, 1, self_idx, self_idx, -(adv_coeff + diff_coeff))
          bc_counter = bc_counter + 1
        else if (ngb_idx == -4) then
          ! Set Dirichlet BCs at top boundary
          call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, "CDS", cps, u, v, ngb_idx)
          call pack_entries(b_coeffs, bc_counter, self_idx, (adv_coeff + diff_coeff)*n_value)
          call pack_entries(mat_coeffs, bc_counter, 1, self_idx, self_idx, -(adv_coeff + diff_coeff))
          bc_counter = bc_counter + 1
        end if
      end do
    end do
    call set_values(b_coeffs, b)
    call set_values(mat_coeffs, M)
    deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)
    deallocate(b_coeffs%idx, b_coeffs%val)
  end subroutine compute_boundary_coeffs

  !subroutine compute_interior_coeffs(mat_coeffs, diagonal_coeffs, boundary_ids, discretisation, u, v)
  !  use types, only: matrix_values
  !  use utils, only: pack_entries

  !  type(matrix_values), intent(inout) :: mat_coeffs
  !  type(matrix_values), intent(inout) :: diagonal_coeffs
  !  integer(accs_int), dimension(:,:), intent(inout) :: boundary_ids
  !  character(len=3), intent(in) :: discretisation
  !  real(accs_real), dimension(:,:), intent(in) :: u, v

  !  integer(accs_int) :: self_idx, ngb_idx
  !  integer(accs_int) :: j
  !  integer(accs_int) :: bc_counter
  !  integer(accs_int) :: mat_assigned
  !  real(accs_real) :: face_area
  !  real(accs_real) :: diff_coeff
  !  real(accs_real) :: adv_coeff

  !  bc_counter = 1
  !  mat_assigned = 1
  !  do self_idx = 1, square_mesh%n
  !    ! Calculate contribution from neighbours
  !    do j = 1, square_mesh%nnb(self_idx)
  !      ngb_idx = square_mesh%nbidx(j, self_idx)
  !      face_area = square_mesh%Af(j, self_idx)
  !      call calc_diffusion_coeff(diff_coeff)
  !      if (ngb_idx > 0) then
  !        call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, discretisation, cps, u, v, 0)
  !        call pack_entries(mat_coeffs, mat_assigned, mat_assigned, self_idx, ngb_idx, adv_coeff + diff_coeff)
  !        print *, 'assigning mat ', mat_assigned, ' to ', self_idx, ngb_idx, &
  !                 mat_coeffs%cglob(mat_assigned), mat_coeffs%rglob(mat_assigned)
  !        call pack_entries(diagonal_coeffs, 1, self_idx, self_idx, self_idx, -(adv_coeff + diff_coeff))
  !        mat_assigned = mat_assigned + 1
  !      else
  !        boundary_ids(bc_counter, 1) = self_idx
  !        boundary_ids(bc_counter, 2) = ngb_idx
  !        bc_counter = bc_counter + 1
  !      end if
  !    end do
  !  end do
  !end subroutine compute_interior_coeffs

  !subroutine compute_boundary_coeffs(mat_coeffs, diagonal_coeffs, boundary_coeffs, &
  !                                   boundary_ids, discretisation, u, v, n_bc_cells)
  !  use types, only: matrix_values, vector_values
  !  use utils, only: pack_entries

  !  type(matrix_values), intent(inout) :: mat_coeffs
  !  type(matrix_values), intent(inout) :: diagonal_coeffs
  !  type(vector_values), intent(inout) :: boundary_coeffs
  !  integer(accs_int), dimension(:,:), intent(in) :: boundary_ids
  !  character(len=3), intent(in) :: discretisation
  !  real(accs_real), dimension(:,:), intent(in) :: u, v
  !  integer(accs_int), intent(in) :: n_bc_cells

  !  integer(accs_int) :: self_idx, ngb_idx
  !  integer(accs_int) :: j
  !  integer(accs_int) :: bc_counter
  !  integer(accs_int) :: row, col
  !  real(accs_real) :: face_area
  !  real(accs_real) :: diff_coeff
  !  real(accs_real) :: adv_coeff
  !  real(accs_real) :: BC_value
  !  real(accs_real) :: n_value, w_value

  !  n_value = 0.0_accs_real
  !  w_value = 1.0_accs_real
  !  bc_counter = 1
  !  call calc_diffusion_coeff(diff_coeff)
  !  do j = 1, n_bc_cells
  !    self_idx = boundary_ids(j,1)
  !    ngb_idx = boundary_ids(j,2)
  !    if (ngb_idx == -1) then
  !      ! Set Dirichlet BCs at left boundary
  !      call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, "CDS", cps, u, v, ngb_idx)
  !      call calc_cell_coords(self_idx, cps, row, col)
  !      BC_value = -(1.0_accs_real - real(row, kind=accs_real)/real(cps, kind=accs_real)) * w_value
  !      call pack_entries(boundary_coeffs, bc_counter, self_idx, (adv_coeff + diff_coeff)*BC_value)
  !      bc_counter = bc_counter + 1
  !      call pack_entries(diagonal_coeffs, 1, self_idx, self_idx, self_idx, -(adv_coeff + diff_coeff))
  !    else if (ngb_idx == -4) then
  !      ! Set Dirichlet BCs at top boundary
  !      call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, "CDS", cps, u, v, ngb_idx)
  !      call pack_entries(boundary_coeffs, bc_counter, self_idx, (adv_coeff + diff_coeff)*n_value)
  !      bc_counter = bc_counter + 1
  !      call pack_entries(diagonal_coeffs, 1, self_idx, self_idx, self_idx, -(adv_coeff + diff_coeff))
  !    end if
  !  end do
  !end subroutine compute_boundary_coeffs

  ! Calculates advection coefficient for neighbouring cell 
  subroutine calc_advection_coeff(ngb_idx, self_idx, face_area, coeff, discretisation, cps, u, v, BC)
    integer(accs_int), intent(in) :: ngb_idx, self_idx
    real(accs_real), intent(in) :: face_area
    real(accs_real), intent(inout) :: coeff
    character(len=3), intent(in) :: discretisation
    integer(accs_int), intent(in) :: cps
    real(accs_real), dimension(:,:) :: u, v
    integer(accs_int), intent(in) :: BC

    integer(accs_int) :: ngb_row, ngb_col       ! neighbour coordinates within grid
    integer(accs_int) :: self_row, self_col     ! cell coordinates within grid

    ! Find where we are in the grid first
    call calc_cell_coords(ngb_idx, cps, ngb_row, ngb_col)
    call calc_cell_coords(self_idx, cps, self_row, self_col)

    coeff = calc_mass_flux(face_area, u, v, ngb_row, ngb_col, self_row, self_col, BC)
    if (discretisation == "UDS") then
      coeff = min(coeff, 0.0_accs_real)
    end if
  end subroutine calc_advection_coeff

  ! Sets diffusion coefficient
  subroutine calc_diffusion_coeff(coeff)
    real(accs_real), intent(inout) :: coeff

    real(accs_real), parameter :: diffusion_factor = 1.e-2

    coeff = -2.*diffusion_factor
  end subroutine calc_diffusion_coeff

  ! Calculates mass flux across given edge. Note: assuming rho = 1 and uniform grid
  function calc_mass_flux(edge_len, u, v, ngb_row, ngb_col, self_row, self_col, BC_flag) result(flux)
    real(accs_real), intent(in) :: edge_len
    real(accs_real), dimension(:,:), intent(in) :: u, v
    integer(accs_int), intent(in) :: ngb_row, ngb_col
    integer(accs_int), intent(in) :: self_row, self_col
    integer(accs_int), intent(in) :: BC_flag

    real(accs_real) :: flux

    flux = 0.

    if (BC_flag == 0) then
      if (ngb_col - self_col == 1) then
        flux = 0.25*(u(ngb_col, ngb_row) + u(self_col, self_row)) * edge_len
      else if (ngb_col - self_col == -1) then
        flux = -0.25*(u(ngb_col, ngb_row) + u(self_col, self_row)) * edge_len
      else if (ngb_row - self_row == 1) then
        flux = 0.25*(v(ngb_col, ngb_row) + v(self_col, self_row)) * edge_len
      else 
        flux = -0.25*(v(ngb_col, ngb_row) + v(self_col, self_row)) * edge_len
      end if
    else if (BC_flag == -1 .or. BC_flag == -2) then
      flux = u(self_col, self_row) * edge_len
    else if (BC_flag == -3 .or. BC_flag == -4) then
      flux = v(self_col, self_row) * edge_len
    else
      print *, 'invalid BC flag'
      stop
    end if
  end function calc_mass_flux

  ! Make the square mesh, and setup ICs for u, v, and our scalar field
  subroutine initialise_scalar_advection(par_env, u, v)

    class(parallel_environment) :: par_env
    real(accs_real), dimension(:,:), allocatable, intent(inout) :: u, v

    integer(accs_int) :: i

    !square_mesh = build_square_mesh(cps, 1.0_accs_real, par_env)
    square_mesh = build_square_mesh(cps, 1.0, par_env)

    ! Allocate velocity and scalar field arrays
    allocate(u(cps,cps))
    allocate(v(cps,cps))

    ! Set IC velocity and scalar fields
    !u = 1.0_accs_real
    !v = 0.0_accs_real
    do i = 1, cps
      u(i,:) = float(i)/float(cps)
      v(:,i) = -float(i)/float(cps)
    end do
  end subroutine initialise_scalar_advection

  ! Assigns source vector
  ! Calculates the row and column indices from flattened vector index
  ! Note: assumes square mesh
  subroutine calc_cell_coords(idx, cps, row, col)
    integer(accs_int), intent(in) :: idx, cps
    integer(accs_int), intent(out) :: row, col

    col = modulo(idx-1,cps) + 1 
    row = (idx-1)/cps + 1
  end subroutine calc_cell_coords

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
