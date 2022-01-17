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
  call compute_fluxes(M, source, u, v)

  call begin_update(M) ! Start the parallel assembly for M

  call begin_update(source) ! Start the parallel assembly for source
  call end_update(M) ! Complete the parallel assembly for M
  call end_update(source) ! Complete the parallel assembly for source

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

  subroutine compute_fluxes(M, b, u, v)
    use constants, only : insert_mode, add_mode
    use types, only : matrix_values, vector_values
    use utils, only : set_values, pack_entries

    class(matrix), intent(inout) :: M   
    class(vector), intent(inout) :: b   
    real(accs_real), dimension(:,:), intent(in) :: u, v
    
    integer(accs_int) :: cps
    integer(accs_int) :: n_bc_cells, n_int_cells
    
    cps = int(sqrt(real(square_mesh%n)))

    ! Loop over cells computing advection and diffusion fluxes
    n_int_cells = calc_matrix_nnz()
    call compute_interior_coeffs(M, "CDS", u, v, n_int_cells)

    ! Loop over boundaries
    n_bc_cells = calc_rhs_nnz(cps)
    call compute_boundary_coeffs(M, b, "CDS", u, v, n_bc_cells)

  end subroutine compute_fluxes

  ! Note: this assumes a 2d grid
  pure function calc_matrix_nnz() result(nnz)
    integer(accs_int) :: nnz

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
          mat_counter = mat_counter + 1
          adv_coeff_total = adv_coeff_total + adv_coeff
          diff_coeff_total = diff_coeff_total + diff_coeff
        else
          call pack_entries(mat_coeffs, 1, mat_counter, self_idx, -1, 0.0_accs_real)
          mat_counter = mat_counter + 1
        end if
      end do
      call pack_entries(mat_coeffs, 1, mat_counter, self_idx, self_idx, -(adv_coeff_total + diff_coeff_total))
      mat_counter = mat_counter + 1
      call set_values(mat_coeffs, mat)
    end do

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
          call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, discretisation, cps, u, v, ngb_idx)
          call calc_cell_coords(self_idx, cps, row, col)
          BC_value = -(1.0_accs_real - real(row, kind=accs_real)/real(cps, kind=accs_real)) * w_value
          call pack_entries(b_coeffs, bc_counter, self_idx, (adv_coeff + diff_coeff)*BC_value)
          call pack_entries(mat_coeffs, bc_counter, 1, self_idx, self_idx, -(adv_coeff + diff_coeff))
          bc_counter = bc_counter + 1
        else if (ngb_idx == -4) then
          ! Set Dirichlet BCs at top boundary
          call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, discretisation, cps, u, v, ngb_idx)
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