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
  integer(accs_int) :: i

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

  !! Evaluate right-hand-side vector
  !call set_zero(source)

  call begin_update(source) ! Start the parallel assembly for source
  call end_update(M) ! Complete the parallel assembly for M
  call end_update(source) ! Complete the parallel assembly for source

  !! Modify matrix and right-hand-side vector to apply Dirichlet boundary conditions
  !call apply_dirichlet_bcs(M, source, -1)
  !call apply_dirichlet_bcs(M, source, -2)
  !call apply_dirichlet_bcs(M, source, -3)
  !call apply_dirichlet_bcs(M, source, -4)
  !call begin_update(source) ! Start the parallel assembly for source
  !call finalise(M)

  call end_update(scalar) ! Complete the parallel assembly for scalar
  !call end_update(source) ! Complete the parallel assembly for source

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

    type(matrix_values) :: mat_coeffs
    type(vector_values) :: b_coeffs
    
    mat_coeffs%mode = add_mode
    b_coeffs%mode = add_mode

    allocate(mat_coeffs%rglob(1))
    allocate(mat_coeffs%cglob(1))
    allocate(mat_coeffs%val(1))
    allocate(b_coeffs%idx(1))
    allocate(b_coeffs%val(1))

    cps = int(sqrt(real(square_mesh%n)))

    ! Loop over cells computing advection and diffusion fluxes
    ! Note, n_cols corresponds to the number of cells in the mesh
    do self_idx = 1, square_mesh%n
      ! Calculate contribution from neighbours
      diff_coeff_total = 0.
      adv_coeff_total = 0.
      do j = 1, square_mesh%nnb(self_idx)
        ngb_idx = square_mesh%nbidx(j, self_idx)
        call calc_diffusion_coeff(diff_coeff, ngb_idx, self_idx)
        ! ALEXEI: only setting matrix elements that correspond to entries within the domain is the same as zero-gradient BCs when using a uniform square mesh.
        if (ngb_idx > 0) then
          call calc_advection_coeff(ngb_idx, self_idx, square_mesh%Af(j,self_idx), adv_coeff, "CDS", cps, u, v, 0)
          call pack_entries(mat_coeffs, 1, 1, self_idx, ngb_idx, adv_coeff + diff_coeff)
          call set_values(mat_coeffs, M)
          diff_coeff_total = diff_coeff_total + diff_coeff
          adv_coeff_total = adv_coeff_total + adv_coeff
        else if (ngb_idx == -1) then
          ! Set Dirichlet BCs at left boundary
          call calc_advection_coeff(ngb_idx, self_idx, square_mesh%Af(j, self_idx), adv_coeff, "CDS", cps, u, v, ngb_idx)
          call calc_cell_coords(self_idx, cps, row, col)
          BC_value = -(1. - float(row)/float(cps)) * w_value
          call pack_entries(b_coeffs, 1, self_idx, (adv_coeff + diff_coeff)*BC_value)
          call set_values(b_coeffs, b)
          diff_coeff_total = diff_coeff_total + diff_coeff
          adv_coeff_total = adv_coeff_total + adv_coeff
        else if (ngb_idx == -4) then
          ! Set Dirichlet BCs at top boundary
          call calc_advection_coeff(ngb_idx, self_idx, square_mesh%Af(j, self_idx), adv_coeff, "CDS", cps, u, v, ngb_idx)
          call pack_entries(b_coeffs, 1, self_idx, (adv_coeff + diff_coeff)*n_value)
          call set_values(b_coeffs, b)
          diff_coeff_total = diff_coeff_total + diff_coeff
          adv_coeff_total = adv_coeff_total + adv_coeff
        end if
      end do
      call pack_entries(mat_coeffs, 1, 1, self_idx, self_idx, -(adv_coeff_total + diff_coeff_total))
      call set_values(mat_coeffs, M)
    end do

    deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)
    deallocate(b_coeffs%idx, b_coeffs%val)
  end subroutine compute_fluxes

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

  ! Calculates diffusion coefficient
  subroutine calc_diffusion_coeff(coeff, ngb_idx, self_idx)
    real(accs_real), intent(inout) :: coeff
    integer(accs_int), intent(in) :: ngb_idx, self_idx

    real(accs_real), parameter :: diffusion_factor = 1.e-2

    ! Because mesh is assumed to be square and uniform, all coefficients (apart from P) 
    ! have the same value
    if (self_idx .ne. ngb_idx) then
      coeff = -2.*diffusion_factor
    else
      coeff = 8*diffusion_factor
    end if
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

    square_mesh = build_square_mesh(cps, 1.0_accs_real, par_env)

    ! Allocate velocity and scalar field arrays
    allocate(u(cps,cps))
    allocate(v(cps,cps))

    ! Set IC velocity and scalar fields
    !v = 1.0_accs_real
    !u = 0.0_accs_real
    do i = 1, cps
      u(i,:) = float(i)/float(cps)
      v(:,i) = -float(i)/float(cps)
    end do
  end subroutine initialise_scalar_advection

  ! Assigns vector to zero
  subroutine set_zero(vec)
    use constants, only : add_mode
    use types, only : vector_values
    use utils, only : set_values, pack_entries

    class(vector), intent(inout) :: vec
    type(vector_values) :: vec_data
    integer(accs_int) :: i, nloc

    vec_data%mode = add_mode
    nloc = square_mesh%n
    allocate(vec_data%idx(nloc))
    allocate(vec_data%val(nloc))

    do i = 1, nloc
      call pack_entries(vec_data, i, i, 0.0_accs_real)
    end do

    call set_values(vec_data, vec)
  end subroutine set_zero

  ! Assigns source vector
  subroutine compute_source_terms(source)
    use constants, only : add_mode
    use types, only : vector_values
    use utils, only : set_values, pack_entries

    class(vector), intent(inout) :: source
    type(vector_values) :: data
    integer(accs_int) :: i, n_non_zero, cps, idx, start

    data%mode = add_mode
    cps = int(sqrt(real(square_mesh%n)))
    n_non_zero = cps
    allocate(data%idx(n_non_zero))
    allocate(data%val(n_non_zero))

    call set_zero(source)

    start = 0
    do i = 1, n_non_zero
      idx = calc_flat_array_index(cps, i+start, cps)
      call pack_entries(data, i, idx, 1.0_accs_real)
    end do
    call set_values(data, source)
  end subroutine compute_source_terms

  ! Sets values of the scalar field on the boundaries
  subroutine set_scalar_BCs(scalar)
    use constants, only : add_mode
    use types, only : vector_values
    use utils, only : set_values, pack_entries

    class(vector), intent(inout) :: scalar
    type(vector_values) :: data
    integer(accs_int) :: i, j, n_non_zero, nz_grid_size, cps, start, idx, data_idx

    data%mode = add_mode
    nz_grid_size = 10
    n_non_zero = nz_grid_size*nz_grid_size
    allocate(data%idx(n_non_zero))
    allocate(data%val(n_non_zero))

    call set_zero(scalar)

    cps = int(sqrt(real(square_mesh%n)))
    start = int(cps/2 - nz_grid_size/2)
    do i = 1, nz_grid_size
      do j = 1, nz_grid_size
        idx = calc_flat_array_index(j+start, i+start, cps)
        data_idx = calc_flat_array_index(j, i, nz_grid_size)
        call pack_entries(data, data_idx, idx, real(square_mesh%n/n_non_zero, kind=accs_real))
      end do
    end do
    call set_values(data, scalar)
  end subroutine set_scalar_BCs

  ! Calculates the row and column indices from flattened vector index
  ! Note: assumes square mesh
  subroutine calc_cell_coords(idx, cps, row, col)
    integer(accs_int), intent(in) :: idx, cps
    integer(accs_int), intent(out) :: row, col

    col = modulo(idx-1,cps) + 1 
    row = (idx-1)/cps + 1
  end subroutine calc_cell_coords

  ! Computes the index of a flat array from the provided row, col and array side length
  ! Note: assumes square array
  pure function calc_flat_array_index(row, col, cps) result(idx)
    integer(accs_int), intent(in) :: row, col, cps
    integer(accs_int) :: idx

    idx = (row-1)*cps + col
  end function calc_flat_array_index

  subroutine write_data(filename, solution)
    character(len=*), intent(in) :: filename
    real(accs_real), dimension(:), intent(in) :: solution

    open(10, file=filename, form='unformatted')
    write(10) solution
    close(10)
  end subroutine write_data
  
  subroutine apply_dirichlet_bcs(M, b, boundary)

    use constants, only : add_mode
    use mat, only : set_eqn
    use types, only : vector_values, matrix_values, matrix, vector, mesh
    use utils, only : set_values, pack_entries
    use kinds, only: accs_int, accs_real
  
    implicit none
  
    class(matrix), intent(inout) :: M
    class(vector), intent(inout) :: b
    integer(accs_int), intent(in) :: boundary
  
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

              if (nbidx == boundary) then
                boundary_coeff = (2.0 / square_mesh%h) * square_mesh%Af(j, i)
                boundary_val = compute_boundary_val(boundary, j)

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

  pure function compute_boundary_val(boundary, i) result(val)
    integer(accs_int), intent(in) :: boundary
    integer(accs_int), intent(in) :: i
    real(accs_real) :: val
    integer(accs_int) :: cps 

    cps = int(sqrt(real(square_mesh%n)))

    if ((boundary == -1 .or. boundary == -3) .and. i > cps*0.25 .and. i < cps*0.75) then
      val = 2.0_accs_real
    end if
  end function compute_boundary_val

!----------------------------------------------------------------------!

  !subroutine eval_rhs(b)

  !  use constants, only : add_mode
  !  use types, only : vector_values
  !  use utils, only : set_values, pack_entries
  !  
  !  class(vector), intent(inout) :: b

  !  integer(accs_int) :: i
  !  integer(accs_int) :: nloc
  !  real(accs_real) :: h
  !  real(accs_real) :: r

  !  type(vector_values) :: val_dat
  !  
  !  val_dat%mode = add_mode
  !  allocate(val_dat%idx(1))
  !  allocate(val_dat%val(1))

  !  nloc = square_mesh%nlocal

  !  h = square_mesh%h

  !  ! this is currently setting 1 vector value at a time
  !  ! consider changing to doing all the updates in one go
  !  ! to do only 1 call to eval_cell_rhs and set_values
  !  associate(idx => square_mesh%idx_global)
  !    do i = 1, nloc
  !      associate(x => square_mesh%xc(1, i), &
  !           y => square_mesh%xc(2, i), &
  !           V => square_mesh%vol(i))
  !        call eval_cell_rhs(x, y, h**2, r)
  !        r = V * r
  !        call pack_entries(val_dat, 1, idx(i), r)
  !        call set_values(val_dat, b)
  !      end associate
  !    end do
  !  end associate
  !  
  !  deallocate(val_dat%idx)
  !  deallocate(val_dat%val)
  !  
  !end subroutine eval_rhs

  !!> @brief Apply forcing function
  !pure subroutine eval_cell_rhs (x, y, H, r)
  !  
  !  real(accs_real), intent(in) :: x, y, H
  !  real(accs_real), intent(out) :: r
  !  
  !  r = 0.0 &
  !       + 0 * (x + y + H) ! Silence unused dummy argument error
  !  
  !end subroutine eval_cell_rhs

  !subroutine discretise_poisson(M)

  !  use constants, only : insert_mode
  !  use types, only : matrix_values
  !  use utils, only : set_values, pack_entries
  !  
  !  class(matrix), intent(inout) :: M

  !  type(matrix_values) :: mat_coeffs
  !  integer(accs_int) :: i, j

  !  integer(accs_int) :: row, col
  !  real(accs_real) :: coeff_f, coeff_p, coeff_nb
  !  
  !  mat_coeffs%mode = insert_mode

  !  !! Loop over cells
  !  do i = 1, square_mesh%nlocal
  !    !> @todo: Doing this in a loop is awful code - malloc maximum coefficients per row once,
  !    !!        filling from front, and pass the number of coefficients to be set, requires
  !    !!        modifying the matrix_values type and the implementation of set_values applied to
  !    !!        matrices.
  !    associate(idxg=>square_mesh%idx_global(i), &
  !              nnb=>square_mesh%nnb(i))
  !      
  !      allocate(mat_coeffs%rglob(1))
  !      allocate(mat_coeffs%cglob(1 + nnb))
  !      allocate(mat_coeffs%val(1 + nnb))

  !      row = idxg
  !      coeff_p = 0.0_accs_real
  !    
  !      !! Loop over faces
  !      do j = 1, nnb
  !        coeff_f = (1.0 / square_mesh%h) * square_mesh%Af(j, i)
  !        associate(nbidxg=>square_mesh%nbidx(j, i))

  !          if (nbidxg > 0) then
  !            !! Interior face
  !            coeff_p = coeff_p - coeff_f
  !            coeff_nb = coeff_f
  !            col = nbidxg
  !          else
  !            col = -1
  !            coeff_nb = 0.0_accs_real
  !          end if
  !          call pack_entries(mat_coeffs, 1, j + 1, row, col, coeff_nb)

  !        end associate
  !      end do

  !      !! Add the diagonal entry
  !      col = row
  !      call pack_entries(mat_coeffs, 1, 1, row, col, coeff_p)

  !      !! Set the values
  !      call set_values(mat_coeffs, M)

  !      deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)
  !      
  !    end associate

  !  end do
  !  
  !end subroutine discretise_poisson

  !subroutine apply_dirichlet_bcs(M, b)

  !  use constants, only : add_mode
  !  use mat, only : set_eqn
  !  use types, only : vector_values, matrix_values, matrix, vector, mesh
  !  use utils, only : set_values, pack_entries
  !  use kinds, only: accs_int, accs_real
  !
  !  implicit none
  !
  !  class(matrix), intent(inout) :: M
  !  class(vector), intent(inout) :: b
  !
  !  integer(accs_int) :: i, j
  !  real(accs_real) :: boundary_coeff, boundary_val

  !  integer(accs_int) :: idx, row, col
  !  real(accs_real) :: r, coeff
  !  
  !  type(vector_values) :: vec_values
  !  type(matrix_values) :: mat_coeffs
  !
  !  allocate(mat_coeffs%rglob(1))
  !  allocate(mat_coeffs%cglob(1))
  !  allocate(mat_coeffs%val(1))
  !  allocate(vec_values%idx(1))
  !  allocate(vec_values%val(1))

  !  mat_coeffs%mode = add_mode
  !  vec_values%mode = add_mode

  !  associate(idx_global=>square_mesh%idx_global)
  !
  !    do i = 1, square_mesh%nlocal
  !      if (minval(square_mesh%nbidx(:, i)) < 0) then
  !        coeff = 0.0_accs_real 
  !        r = 0.0_accs_real
  !        
  !        row = idx_global(i)
  !        col = idx_global(i)
  !        idx = idx_global(i)

  !        do j = 1, square_mesh%nnb(i)

  !          associate(nbidx=>square_mesh%nbidx(j, i))

  !            if (nbidx < 0) then
  !              boundary_coeff = (2.0 / square_mesh%h) * square_mesh%Af(j, i)
  !              boundary_val = rhs_val(i, j)

  !              ! Coefficient
  !              coeff = coeff - boundary_coeff

  !              ! RHS vector
  !              r = r - boundary_val * boundary_coeff
  !            end if

  !          end associate
  !        end do

  !        call pack_entries(mat_coeffs, 1, 1, row, col, coeff)
  !        call pack_entries(vec_values, 1, idx, r)

  !        call set_values(mat_coeffs, M)
  !        call set_values(vec_values, b)
  !        
  !      end if
  !    end do
  !
  !  end associate
  !
  !  deallocate(vec_values%idx)
  !  deallocate(vec_values%val)
  !
  !end subroutine apply_dirichlet_bcs

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

  !subroutine initialise_poisson(par_env)

  !  class(parallel_environment) :: par_env

  !  square_mesh = build_square_mesh(cps, 1.0_accs_real, par_env)
  !  
  !end subroutine initialise_poisson

  !pure function rhs_val(i, f) result(r)

  !  integer(accs_int), intent(in) :: i !> Cell index
  !  integer(accs_int), intent(in), optional :: f !> Face index (local wrt cell)

  !  real(accs_real) :: r

  !  if (present(f)) then
  !    !! Face-centred value
  !    associate(y => square_mesh%xf(2, f, i))
  !      r = y
  !    end associate
  !  else
  !    !! Cell-centred value
  !    associate(y => square_mesh%xc(2, i))
  !      r = y
  !    end associate
  !  end if
  !  
  !end function rhs_val

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
