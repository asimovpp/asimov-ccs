!> @brief Test that the flux matrix has been computed correctly
!
!> @description Compares the matrix calculated for flows in the +x and +y directions with
!> central and upwind differencing to the known matrix
program test_compute_fluxes

  use testing_lib
  use types, only: field, upwind_field, central_field
  use mesh_utils, only : build_square_mesh
  use fv, only: compute_fluxes, calc_cell_coords
  use utils, only : update, initialise, &
                set_global_size, pack_entries, set_values
  use vec, only : create_vector
  use mat, only : create_matrix, set_nnz
  use solver, only : axpy, norm
  use constants, only: add_mode, insert_mode
  use bc_constants

  implicit none

  type(mesh) :: square_mesh
  type(bc_config) :: bcs
  type(vector_init_data) :: vec_sizes
  class(field), allocatable :: scalar
  class(field), allocatable :: u, v
  integer(accs_int), parameter :: cps = 5
  integer(accs_int) :: direction, discretisation
  integer, parameter :: x_dir = 1, y_dir = 2
  integer, parameter :: central = -1, upwind = -2

  call init()

  square_mesh = build_square_mesh(cps, 1.0_accs_real, par_env)

  bcs%region(1) = bc_region_left
  bcs%region(2) = bc_region_right
  bcs%region(3) = bc_region_top
  bcs%region(4) = bc_region_bottom
  bcs%bc_type(:) = bc_type_dirichlet
  bcs%endpoints(:,:) = 1.0_accs_real
    
  do direction = x_dir, y_dir
    do discretisation = upwind, central
    
      if (discretisation == central) then
        allocate(central_field :: scalar)
        allocate(central_field :: u)
        allocate(central_field :: v)
      else if (discretisation == upwind) then
        allocate(upwind_field :: scalar)
        allocate(upwind_field :: u)
        allocate(upwind_field :: v)
      else
        write(message, *) 'Invalid discretisation type selected'
        call stop_test(message)
      end if

      call initialise(vec_sizes)
      call set_global_size(vec_sizes, square_mesh, par_env)
      call create_vector(vec_sizes, scalar%vec)
      call create_vector(vec_sizes, u%vec)
      call create_vector(vec_sizes, v%vec)

      call set_velocity_fields(square_mesh, direction, u, v)
      call run_compute_fluxes_test(scalar, u, v, bcs, square_mesh, cps, direction, discretisation)
      call tidy_velocity_fields(scalar, u, v)
    end do
  end do

  call fin()

  contains

  !> @brief Sets the velocity field in the desired direction and discretisation
  !
  !> @param[in] cell_mesh - The mesh structure
  !> @param[in] direction - Integer indicating the direction of the velocity field
  !> @param[out] u, v     - The velocity fields in x and y directions
  subroutine set_velocity_fields(cell_mesh, direction, u, v)
    use meshing, only: set_cell_location, get_global_index
    class(mesh), intent(in) :: cell_mesh
    integer(accs_int), intent(in) :: direction
    class(field), intent(inout), allocatable :: u, v
    type(cell_locator) :: self_loc
    type(vector_values) :: u_vals, v_vals
    integer(accs_int) :: local_idx, self_idx
    real(accs_real) :: u_val, v_val

    u_vals%mode = insert_mode
    v_vals%mode = insert_mode
    
    associate(n_local => cell_mesh%nlocal)
      allocate(u_vals%idx(n_local))
      allocate(v_vals%idx(n_local))
      allocate(u_vals%val(n_local))
      allocate(v_vals%val(n_local))
      
      ! Set IC velocity fields
      do local_idx = 1, n_local
        call set_cell_location(self_loc, cell_mesh, local_idx)
        call get_global_index(self_loc, self_idx)

        if (direction == x_dir) then
          u_val = 1.0_accs_real
          v_val = 0.0_accs_real
        else if (direction == y_dir) then
          u_val = 0.0_accs_real
          v_val = 1.0_accs_real
        end if

        call pack_entries(u_vals, local_idx, self_idx, u_val)
        call pack_entries(v_vals, local_idx, self_idx, v_val)
      end do
    end associate
    call set_values(u_vals, u%vec)
    call set_values(v_vals, v%vec)

    deallocate(u_vals%idx, v_vals%idx, u_vals%val, v_vals%val)
  end subroutine set_velocity_fields

  !> @brief Deallocates the velocity fields
  !
  !> @param[in] scalar - The scalar field structure
  !> @param[in] u, v   - The velocity fields to deallocate
  subroutine tidy_velocity_fields(scalar, u, v)
    class(field), allocatable :: scalar
    class(field), allocatable :: u, v

    deallocate(scalar)
    deallocate(u)
    deallocate(v)
  end subroutine tidy_velocity_fields

  !> @brief Compares the matrix computed for a given velocity field and discretisation to the known solution
  !
  !> @param[in] scalar         - The scalar field structure
  !> @param[in] u, v           - The velocity field structures
  !> @param[in] bcs            - The BC structure
  !> @param[in] cell_mesh      - The mesh structure
  !> @param[in] cps            - The number of cells per side in the (square) mesh 
  !> @param[in] flow_direction - Integer indicating the direction of the flow 
  !> @param[in] discretisation - Integer indicating the discretisation scheme being tested 
  subroutine run_compute_fluxes_test(scalar, u, v, bcs, cell_mesh, cps, flow_direction, discretisation)
    class(field), intent(in) :: scalar
    class(field), intent(in) :: u, v
    class(bc_config), intent(in) :: bcs
    type(mesh), intent(in) :: cell_mesh
    integer(accs_int), intent(in) :: cps
    integer(accs_int), intent(in) :: flow_direction
    integer(accs_int), intent(in) :: discretisation

    class(matrix), allocatable :: M, M_exact
    class(vector), allocatable :: b, b_exact
    type(vector_init_data) :: vec_sizes
    type(matrix_init_data) :: mat_sizes
    real(accs_real) :: error
    
    call initialise(mat_sizes)
    call initialise(vec_sizes)
    call set_global_size(mat_sizes, cell_mesh, par_env)
    call set_global_size(vec_sizes, cell_mesh, par_env)
    call set_nnz(mat_sizes, 5)
    call create_matrix(mat_sizes, M)
    call create_vector(vec_sizes, b)
    call create_matrix(mat_sizes, M_exact)
    call create_vector(vec_sizes, b_exact)

    call compute_fluxes(scalar, u, v, cell_mesh, bcs, cps, M, b)

    call update(M)
    call update(b)

    call compute_exact_matrix(cell_mesh, flow_direction, discretisation, cps, M_exact, b_exact)

    call update(M_exact)
    call update(b_exact)

    call axpy(-1.0_accs_real, M_exact, M)
    error = norm(M, 1)

    if (error .ge. eps) then
      write(message, *) 'FAIL: matrix difference norm too large ', error
      call stop_test(message)
    end if
    
    call axpy(-1.0_accs_real, b_exact, b)
    error = norm(b, 2)

    if (error .ge. eps) then
      write(message, *) 'FAIL: vector difference norm too large ', error
      call stop_test(message)
    end if

    deallocate(M)
    deallocate(b)
    deallocate(M_exact)
    deallocate(b_exact)
  end subroutine run_compute_fluxes_test

  !> @brief Computes the known flux matrix for the given flow and discretisation
  !
  !> @param[in] cell_mesh      - The (square) mesh
  !> @param[in] flow           - Integer indicating flow direction
  !> @param[in] discretisation - Integer indicating the discretisation scheme being used
  !> @param[in] cps            - Number of cells per side in mesh
  !> @param[out] M             - The resulting matrix
  !> @param[out] b             - The resulting RHS vector
  subroutine compute_exact_matrix(cell_mesh, flow, discretisation, cps, M, b)
    class(mesh), intent(in) :: cell_mesh
    integer(accs_int), intent(in) :: flow
    integer(accs_int), intent(in) :: discretisation
    integer(accs_int), intent(in) :: cps
    class(matrix), allocatable :: M
    class(vector), allocatable :: b

    type(matrix_init_data) :: mat_sizes
    type(matrix_values) :: mat_coeffs
    type(vector_init_data) :: vec_sizes
    type(vector_values) :: vec_coeffs
    real(accs_real) :: diff_coeff, adv_coeff
    real(accs_real) :: proc_zero_factor ! set to zero for non-zero ranks so that we're not double counting when running in parallel
    integer(accs_int) :: i, j
    integer(accs_int) :: row, col
    integer(accs_int) :: mat_counter
    integer(accs_int) :: vec_counter

    call initialise(mat_sizes)
    call set_global_size(mat_sizes, cell_mesh, par_env)
    call set_nnz(mat_sizes, 5)
    call create_matrix(mat_sizes, M)
    
    call initialise(vec_sizes)
    call set_global_size(vec_sizes, cell_mesh, par_env)
    call create_vector(vec_sizes, b)

    mat_coeffs%mode = add_mode
    vec_coeffs%mode = add_mode
    
    allocate(mat_coeffs%rglob(1))
    allocate(mat_coeffs%cglob(2))
    allocate(mat_coeffs%val(2))

    if (par_env%proc_id == 0) then  ! We only want to assign the values once, and there's no point parallelising the test
      proc_zero_factor = 1.0_accs_real
    else 
      proc_zero_factor = 0.0_accs_real
    end if

    ! Advection coefficients first
    if (flow == x_dir .and. discretisation == central) then
      ! CDS and flow along +x direction
      do i = 1, cell_mesh%nglobal
        mat_counter = 1
        if (mod(i, cps) == 1) then
          call pack_entries(mat_coeffs, 1, mat_counter, i, i, -0.3_accs_real*proc_zero_factor) ! Make this more flexible so that the coeffs depend on cps
          mat_counter = mat_counter + 1
          call pack_entries(mat_coeffs, 1, mat_counter, i, i+1, 0.1_accs_real*proc_zero_factor)
          mat_counter = mat_counter + 1
        else if (mod(i, cps) == 0) then
          call pack_entries(mat_coeffs, 1, mat_counter, i, i, -0.1_accs_real*proc_zero_factor)
          mat_counter = mat_counter + 1
          call pack_entries(mat_coeffs, 1, mat_counter, i, i-1, -0.1_accs_real*proc_zero_factor)
          mat_counter = mat_counter + 1
        else
          call pack_entries(mat_coeffs, 1, mat_counter, i, i+1, 0.1_accs_real*proc_zero_factor)
          mat_counter = mat_counter + 1
          call pack_entries(mat_coeffs, 1, mat_counter, i, i-1, -0.1_accs_real*proc_zero_factor)
          mat_counter = mat_counter + 1
        end if
        call set_values(mat_coeffs, M)
      end do
    else if (flow == y_dir .and. discretisation == central) then
      ! CDS and flow along +y direction
      do i = 1, cell_mesh%nglobal
        mat_counter = 1
        if (i .le. cps) then
          call pack_entries(mat_coeffs, 1, mat_counter, i, i, -0.3_accs_real*proc_zero_factor)
          mat_counter = mat_counter + 1
        else if (i > cell_mesh%nglobal - cps) then
          call pack_entries(mat_coeffs, 1, mat_counter, i, i, -0.1_accs_real*proc_zero_factor)
          mat_counter = mat_counter + 1
        end if

        if (i + cps .le. cell_mesh%nglobal) then
          call pack_entries(mat_coeffs, 1, mat_counter, i, i+cps, 0.1_accs_real*proc_zero_factor)
          mat_counter = mat_counter + 1
        end if
        if (i - cps > 0) then
          call pack_entries(mat_coeffs, 1, mat_counter, i, i-cps, -0.1_accs_real*proc_zero_factor)
          mat_counter = mat_counter + 1
        end if
        call set_values(mat_coeffs, M)
      end do
    else if (flow == x_dir .and. discretisation == upwind) then
      ! UDS and flow along +x direction
      do i = 1, cell_mesh%nglobal
        mat_counter = 1
        if (mod(i, cps) .ne. 1) then
          call pack_entries(mat_coeffs, 1, mat_counter, i, i, 0.2_accs_real*proc_zero_factor)
          mat_counter = mat_counter + 1
          call pack_entries(mat_coeffs, 1, mat_counter, i, i-1, -0.2_accs_real*proc_zero_factor)
          mat_counter = mat_counter + 1
          call set_values(mat_coeffs, M)
        end if
      end do
    else if (flow == y_dir .and. discretisation == upwind) then
      ! UDS and flow along +y direction
      do i = cps+1, cell_mesh%nglobal
        mat_counter = 1
        call pack_entries(mat_coeffs, 1, mat_counter, i, i, 0.2_accs_real*proc_zero_factor)
        mat_counter = mat_counter + 1
        call pack_entries(mat_coeffs, 1, mat_counter, i, i-cps, -0.2_accs_real*proc_zero_factor)
        mat_counter = mat_counter + 1
        call set_values(mat_coeffs, M)
      end do
    end if

    deallocate(mat_coeffs%cglob)
    deallocate(mat_coeffs%val)

    allocate(mat_coeffs%cglob(5))
    allocate(mat_coeffs%val(5))

    diff_coeff = -0.01_accs_real
    ! Diffusion coefficients
    do i = 1, cell_mesh%nglobal
      mat_counter = 1
      call pack_entries(mat_coeffs, 1, mat_counter, i, i, -4*diff_coeff*proc_zero_factor)
      mat_counter = mat_counter + 1

      if (i - 1 > 0 .and. mod(i, cps) .ne. 1) then
        call pack_entries(mat_coeffs, 1, mat_counter, i, i-1, diff_coeff*proc_zero_factor)
        mat_counter = mat_counter + 1
      end if
      if (i - cps > 0) then
        call pack_entries(mat_coeffs, 1, mat_counter, i, i-cps, diff_coeff*proc_zero_factor)
        mat_counter = mat_counter + 1
      end if

      if (i + 1 .le. cell_mesh%nglobal .and. mod(i, cps) .ne. 0) then
        call pack_entries(mat_coeffs, 1, mat_counter, i, i+1, diff_coeff*proc_zero_factor)
        mat_counter = mat_counter + 1
      end if
      if (i + cps .le. cell_mesh%nglobal) then
        call pack_entries(mat_coeffs, 1, mat_counter, i, i+cps, diff_coeff*proc_zero_factor)
        mat_counter = mat_counter + 1
      end if

      if (mat_counter < 6) then
        do j = mat_counter, cps
          call pack_entries(mat_coeffs, 1, mat_counter, i, -1, 0.0_accs_real)
          mat_counter = mat_counter + 1
        end do
      end if
      call set_values(mat_coeffs, M)
    end do

    deallocate(mat_coeffs%rglob)
    deallocate(mat_coeffs%cglob)
    deallocate(mat_coeffs%val)
    
    ! Now do the RHS
    ! Advection first
    allocate(vec_coeffs%idx(2*cell_mesh%nglobal/cps))
    allocate(vec_coeffs%val(2*cell_mesh%nglobal/cps))

    vec_counter = 1
    if (discretisation == central) then
      adv_coeff = -1.0_accs_real/cps
    else
      adv_coeff = 0.0_accs_real
    endif 

    if (flow == x_dir) then
      do i = 1, cps
        call pack_entries(vec_coeffs, vec_counter, (i-1)*cps + 1, adv_coeff*proc_zero_factor) 
        vec_counter = vec_counter + 1
        call pack_entries(vec_coeffs, vec_counter, i*cps, adv_coeff*proc_zero_factor) 
        vec_counter = vec_counter + 1
      end do
    else
      do i = 1, cps
        call pack_entries(vec_coeffs, vec_counter, i, adv_coeff*proc_zero_factor) 
        vec_counter = vec_counter + 1
        call pack_entries(vec_coeffs, vec_counter, cell_mesh%nglobal - i + 1, adv_coeff*proc_zero_factor) 
        vec_counter = vec_counter + 1
      end do
    end if
    call set_values(vec_coeffs, b)
    
    deallocate(vec_coeffs%idx)
    deallocate(vec_coeffs%val)

    ! And now diffusion
    allocate(vec_coeffs%idx(4*cps))
    allocate(vec_coeffs%val(4*cps))

    vec_counter = 1
    diff_coeff = 0.01_accs_real
    do i = 1, cell_mesh%nglobal
      call calc_cell_coords(i, cps, row, col)
      if (row == 1 .or. row == cps) then
        call pack_entries(vec_coeffs, vec_counter, i, diff_coeff*proc_zero_factor) 
        vec_counter = vec_counter + 1
      end if
      if (col == 1 .or. col == cps) then
        call pack_entries(vec_coeffs, vec_counter, i, diff_coeff*proc_zero_factor) 
        vec_counter = vec_counter + 1
      end if
    end do
    call set_values(vec_coeffs, b)
    
    deallocate(vec_coeffs%idx)
    deallocate(vec_coeffs%val)
    
  end subroutine compute_exact_matrix

end program test_compute_fluxes
