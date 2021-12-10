!> @brief Program file for scalar advection case
!
!

program scalar_advection

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

  integer(accs_int) :: cps = 10 ! Default value for cells per side

  !real(accs_real) :: err_norm

  double precision :: start_time
  double precision :: end_time

  call initialise_parallel_environment(par_env) 
  call read_command_line_arguments()
  call timer(start_time)

  !call initialise_poisson(par_env)

  ! Init ICs (velocities, BC scalar, mesh, etc)
  call initialise_scalar_advection(par_env, u, v)

  !! Initialise with default values
  call initialise(mat_sizes)
  call initialise(vec_sizes)
  call initialise(scalar_linear_system)

  !! Create stiffness matrix
  call set_global_size(mat_sizes, square_mesh%n, square_mesh%n, par_env)
  call set_nnz(mat_sizes, 5) ! ALEXEI: do we also need to have 5 non-zeros here?
  call create_matrix(mat_sizes, M)

  ! Actually compute the values to fill the matrix
  call compute_fluxes(M, u, v)

  call begin_update(M) ! Start the parallel assembly for M

  !! Create right-hand-side and solution vectors
  call set_global_size(vec_sizes, square_mesh%n, par_env)
  call create_vector(vec_sizes, source)
  call create_vector(vec_sizes, solution)
  call create_vector(vec_sizes, scalar)

  call begin_update(scalar) ! Start the parallel assembly for scalar

  ! Update scalar BCs
  call set_scalar_BCs(scalar)

  !! Evaluate right-hand-side vector
  call compute_source_terms(source)

  call begin_update(source) ! Start the parallel assembly for source
  call end_update(M) ! Complete the parallel assembly for M
  call end_update(source) ! Complete the parallel assembly for source

  !! Modify matrix and right-hand-side vector to apply Dirichlet boundary conditions
  !call apply_dirichlet_bcs(M, source)
  !call begin_update(source) ! Start the parallel assembly for source
  !call finalise(M)

  call end_update(scalar) ! Complete the parallel assembly for scalar
  !call end_update(source) ! Complete the parallel assembly for source

  !! Create linear solver & set options
  call set_linear_system(scalar_linear_system, source, scalar, M, par_env)
  call create_solver(scalar_linear_system, scalar_solver)
  call solve(scalar_solver)

  !! Check solution
  ! ALEXEI: Need to know what the analytic solution for this problem is first.
  !call set_exact_sol(ustar)
  !call axpy(-1.0_accs_real, ustar, scalar)

  !err_norm = norm(scalar, 2) * square_mesh%h
  !if (par_env%proc_id == par_env%root) then
  !   print *, "Norm of error = ", err_norm
  !end if
  
  !! Clean up
  deallocate(scalar)
  deallocate(source)
  !deallocate(ustar)
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

  subroutine compute_fluxes(M, u, v)
    use constants, only : insert_mode
    use types, only : matrix_values
    use utils, only : set_values, pack_entries

    class(matrix), intent(inout) :: M   
    real(accs_real), dimension(:,:), intent(in) :: u, v
    
    integer(accs_int) :: i, j           ! Counters
    integer(accs_int) :: n_rows, n_cols 
    integer(accs_int) :: row, col 
    integer(accs_int), parameter :: n_coeffs = 5
    real(accs_real), dimension(n_coeffs) :: adv_coeffs, diff_coeffs

    type(matrix_values) :: mat_coeffs

    mat_coeffs%mode = insert_mode

    allocate(mat_coeffs%rglob(1))
    allocate(mat_coeffs%cglob(1))
    allocate(mat_coeffs%val(1))

    select type (M)
      type is (matrix_petsc) ! ALEXEI: make this independent of solver implementation
      n_rows = int(sqrt(real(square_mesh%n)))
      n_cols = n_rows
      
      ! calculate diffusion coefficients outside the loop because they shouldn't depend on position
      call calc_diffusion_coeffs(diff_coeffs)

      ! Loop over cells computing advection and diffusion fluxes
      ! Note, n_cols corresponds to the number of cells in the mesh
      do i = 1, n_cols
        call calc_advection_coeffs(i, adv_coeffs, "CDS", n_cols, u, v) ! change string comparison to int comparison.

        ! Assign fluxes to matrix. Coefficients have same ordering as eq 4.51 of Ferziger.
        ! This loop finds the locations of the corresponding cell in the matrix and updates 
        ! its value
        do j = 1, n_coeffs
          call get_matrix_ngbs(i, j, n_cols, row, col) ! Should eventually use mesh to find neighbours
          ! Replace this hack setting row and column values to -1 with proper use of neighbours in mesh
          if (row > 0 .and. col > 0) then
            call pack_entries(mat_coeffs, 1, 1, row, col, adv_coeffs(j) + diff_coeffs(j)) ! I think this is only suitable for serial with
                                                                                          ! this choice of indices
            call set_values(mat_coeffs, M)
          end if
        end do
      end do

      class default
        write(*,*) "Unsupported matrix type"
        stop
    end select
    deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)
  end subroutine compute_fluxes

  ! Calculates advection coefficients and stores them in coefficient array using ordering i
  ! of eq 4.51 of Ferziger (i.e. W, S, P, N, E)
  subroutine calc_advection_coeffs(idx, coeffs, discretisation, cps, u, v)
    integer(accs_int), intent(in) :: idx
    real(accs_real), dimension(5), intent(inout) :: coeffs
    character(len=3), intent(in) :: discretisation
    integer(accs_int), intent(in) :: cps
    real(accs_real), dimension(:,:) :: u, v

    real(accs_real) :: edge_len         ! cell edge length
    integer(accs_int) :: row, col       ! coordinates within grid

    ! Find where we are in the grid first
    call calc_cell_coords(idx, cps, row, col)

    edge_len = 1./real(cps)

    select case(discretisation)
      case("UDS")
        coeffs(1) = min(calc_mass_flux("w", edge_len, u, v, row, col), 0.0_accs_real)
        coeffs(2) = min(calc_mass_flux("s", edge_len, u, v, row, col), 0.0_accs_real)
        coeffs(4) = min(calc_mass_flux("n", edge_len, u, v, row, col), 0.0_accs_real)
        coeffs(5) = min(calc_mass_flux("e", edge_len, u, v, row, col), 0.0_accs_real)
        coeffs(3) = -(coeffs(1) + coeffs(2) + coeffs(4) + coeffs(5))
      case("CDS")
        ! NOTE: Assumes uniform grid! (for now)
        coeffs(1) = calc_mass_flux("w", edge_len, u, v, row, col)
        coeffs(2) = calc_mass_flux("s", edge_len, u, v, row, col)
        coeffs(4) = calc_mass_flux("n", edge_len, u, v, row, col)
        coeffs(5) = calc_mass_flux("e", edge_len, u, v, row, col)
        coeffs(3) = -(coeffs(1) + coeffs(2) + coeffs(4) + coeffs(5))
      case default
        write(*,*) "Invalid discretisation scheme provided. Aborting"
        stop
    end select

  end subroutine calc_advection_coeffs

  ! Calculates diffusion coefficients and stores them in coefficient array using ordering i
  ! of eq 4.51 of Ferziger (i.e. W, S, P, N, E)
  subroutine calc_diffusion_coeffs(coeffs)
    real(accs_real), dimension(5), intent(inout) :: coeffs

    real(accs_real), parameter :: Gamma = 0.001

    ! Because mesh is assumed to be square and uniform, all coefficients (apart from P) 
    ! have the same value
    coeffs = -2.*Gamma 
    coeffs(3) = -4*coeffs(1)
  end subroutine calc_diffusion_coeffs

  ! Calculates mass flux across given edge. Note: assuming rho = 1
  function calc_mass_flux(edge, edge_len, u, v, row, col) result(flux)
    character, intent(in) :: edge
    real(accs_real), intent(in) :: edge_len
    real(accs_real), dimension(:,:), intent(in) :: u, v
    integer(accs_int), intent(in) :: row, col
    integer(accs_int) :: n_cols

    real(accs_real) :: flux

    n_cols = int(sqrt(real(square_mesh%n)))

    select case(edge)
      case("w")
        if (col > 1) then
          flux = 0.5*(u(col, row) + u(col-1, row)) * edge_len
        else
          flux = u(col, row) * edge_len
        end if
      case("e")
        if (col < n_cols) then
          flux = 0.5*(u(col, row) + u(col+1, row)) * edge_len
        else
          flux = u(col, row) * edge_len
        end if
      case("s")
        if (row > 1) then
          flux = 0.5*(v(col, row) + v(col, row-1)) * edge_len
        else
          flux = v(col, row) * edge_len
        end if
      case("n")
        if (row < n_cols) then
          flux = 0.5*(v(col, row) + v(col, row+1)) * edge_len
        else
          flux = v(col, row) * edge_len
        end if
      case default
        write(*,*) "Invalid edge provided. Aborting"
        stop
    end select
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
    do i = 1, cps
      u(:,i) = -real(i)/real(cps)
      v(i,:) = real(i)/real(cps)
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
      !vec_data%idx(i) = i
      !vec_data%val(i) = 0
      call pack_entries(vec_data, i, i, 0.0_accs_real)
    end do

    call set_values(vec_data, vec)
  end subroutine set_zero

  ! Assigns source vector
  subroutine compute_source_terms(source)
    class(vector), intent(inout) :: source

    call set_zero(source)
  end subroutine compute_source_terms

  ! Sets values of the scalar field on the boundaries
  subroutine set_scalar_BCs(scalar)
    use constants, only : add_mode
    use types, only : vector_values
    use utils, only : set_values, pack_entries

    class(vector), intent(inout) :: scalar
    type(vector_values) :: data
    integer(accs_int) :: i, cps

    data%mode = add_mode
    cps = int(sqrt(real(square_mesh%n)))
    allocate(data%idx(cps))
    allocate(data%val(cps))

    call set_zero(scalar)

    do i = 1, cps
      !data%idx = i + square_mesh%n - cps
      !data%val = 1
      call pack_entries(data, i, i + square_mesh%n - cps, 1.0_accs_real)
      !print *, 'set_scalar idx ', data%idx
    end do
    call set_values(data, scalar)
  end subroutine set_scalar_BCs

  ! Calculates the row and column indices from flattened vector index
  ! Note: assumes square mesh
  subroutine calc_cell_coords(idx, cps, row, col)
    integer(accs_int), intent(in) :: idx, cps
    integer(accs_int), intent(inout) :: row, col

    col = modulo(idx-1,cps) + 1 
    row = (idx-1)/cps + 1
  end subroutine calc_cell_coords

  ! Calculates the index of a flattened 2d array
  !pure function calc_flat_index(row, col, cps) result(index)
  !  integer(accs_int), intent(in) :: row, col, cps
  !  integer(accs_int) :: index

  !  index = col + (row - 1)*cps
  !end function calc_flat_index

  ! Computes the row and column indices of 
  subroutine get_matrix_ngbs(cell_index, coeff_index, cps, row, col)
    integer(accs_int), intent(in) :: cell_index
    integer(accs_int), intent(in) :: coeff_index
    integer(accs_int), intent(in) :: cps
    integer(accs_int), intent(inout) :: row, col
    integer(accs_int) :: central_col    ! column of central cell
    integer(accs_int) :: n_cols

    n_cols = int(sqrt(real(square_mesh%n)))

    call calc_cell_coords(cell_index, cps, row, central_col)
    if (coeff_index == 3) then
      col = central_col
    else if (coeff_index == 1) then
      if (cell_index-1 < 1) then
        row = -1
        col = -1
      else
        call calc_cell_coords(cell_index-1, cps, row, col)
        col = central_col
      end if
    else if (coeff_index == 2) then
      if (cell_index-cps < 1) then
        row = -1
        col = -1
      else
        call calc_cell_coords(cell_index-cps, cps, row, col)
        col = central_col
      end if
    else if (coeff_index == 4) then
      if (cell_index+cps > n_cols) then
        row = -1
        col = -1
      else
        call calc_cell_coords(cell_index+cps, cps, row, col)
        col = central_col
      end if
    else if (coeff_index == 5) then
      if (cell_index+1 > n_cols) then
        row = -1
        col = -1
      else
        call calc_cell_coords(cell_index+1, cps, row, col)
        col = central_col
      end if
    end if

    if (row > n_cols .or. row < -1 .or. row == 0) then
      print *, 'invalid row ', row
      stop
    end if
    if (col > n_cols .or. col < -1 .or. col == 0) then
      print *, 'invalid col ', col
      stop
    end if

  end subroutine

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

  !subroutine set_exact_sol(ustar)

  !  use constants, only : insert_mode
  !  use types, only : vector_values
  !  use utils, only : set_values, pack_entries
  !  
  !  class(vector), intent(inout) :: ustar

  !  type(vector_values) :: vec_values
  !  integer(accs_int) :: i
  !  
  !  allocate(vec_values%idx(1))
  !  allocate(vec_values%val(1))
  !  vec_values%mode = insert_mode
  !  do i = 1, square_mesh%nlocal
  !     associate(idx=>square_mesh%idx_global(i))
  !       call pack_entries(vec_values, 1, idx, rhs_val(i))
  !       call set_values(vec_values, ustar)
  !     end associate
  !  end do
  !  deallocate(vec_values%idx)
  !  deallocate(vec_values%val)

  !  call update(ustar)
  !end subroutine set_exact_sol

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
