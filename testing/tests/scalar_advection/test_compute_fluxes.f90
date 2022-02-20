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
                set_global_size, pack_entries, set_values, axpy, norm
  use vec, only : create_vector
  use mat, only : create_matrix, set_nnz
  use BC_constants


  type(mesh) :: square_mesh
  type(BC_config) :: BCs
  class(field), allocatable :: u, v
  integer(accs_int), parameter :: cps = 5
  integer(accs_int) :: direction, discretisation

  call init()

  square_mesh = build_square_mesh(cps, 1.0_accs_real, par_env)

  BCs%region(1) = BC_region_left
  BCs%region(2) = BC_region_right
  BCs%region(3) = BC_region_top
  BCs%region(4) = BC_region_bottom
  BCs%BC_type(:) = BC_type_dirichlet
  BCs%endpoints(:,:) = 1.0_accs_real

  do direction = 1, 2
    do discretisation = 1, 2
      call set_velocity_fields(direction, cps, discretisation, u, v)
      call run_compute_fluxes_test(u, v, BCs, square_mesh, cps, direction, discretisation)
      call tidy_velocity_fields(u, v)
    end do
  end do

  call fin()

  contains

  !> @brief Sets the velocity field in the desired direction and discretisation
  !
  !> @param[in] direction      - Integer indicating the direction of the velocity field
  !> @param[in] cps            - Number of cells per side in (square) mesh
  !> @param[in] discretisation - Integer indicating which discretisation scheme to use
  !> @param[out] u, v          - The velocity fields in x and y directions
  subroutine set_velocity_fields(direction, cps, discretisation, u, v)
    integer(accs_int), intent(in) :: direction
    integer(accs_int), intent(in) :: cps
    integer(accs_int), intent(in) :: discretisation
    class(field), intent(inout), allocatable :: u, v

    if (discretisation == 1) then
      allocate(central_field :: u)
      allocate(central_field :: v)
    else if (discretisation == 2) then
      allocate(upwind_field :: u)
      allocate(upwind_field :: v)
    else
      write(message, *) 'Invalid discretisation type selected'
      call stop_test(message)
    end if
    allocate(u%val(cps,cps))
    allocate(v%val(cps,cps))

    if (direction == 1) then
       u%val(:,:) = 1.0_accs_real
       v%val(:,:) = 0.0_accs_real
    else if (direction == 2) then
      u%val(:,:) = 0.0_accs_real
      v%val(:,:) = 1.0_accs_real
    end if
  end subroutine set_velocity_fields

  !> @brief Deallocates the velocity fields
  !
  !> @param[in] u, v - The velocity fields to deallocate
  subroutine tidy_velocity_fields(u, v)
    class(field), allocatable :: u, v

    deallocate(u%val)
    deallocate(v%val)
    deallocate(u)
    deallocate(v)
  end subroutine tidy_velocity_fields

  !> @brief Compares the matrix computed for a given velocity field and discretisation to the known solution
  !
  !> @param[in] u, v           - The velocity fields
  !> @param[in] BCs            - The BC structure
  !> @param[in] cell_mesh      - The mesh structure
  !> @param[in] cps            - The number of cells per side in the (square) mesh 
  !> @param[in] flow_direction - Integer indicating the direction of the flow 
  !> @param[in] discretisation - Integer indicating the discretisation scheme being tested 
  subroutine run_compute_fluxes_test(u, v, BCs, cell_mesh, cps, flow_direction, discretisation)
    class(field), intent(in) :: u, v
    class(BC_config), intent(in) :: BCs
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
    call set_global_size(mat_sizes, cell_mesh%nglobal, cell_mesh%nglobal, par_env)
    call set_global_size(vec_sizes, cell_mesh%nglobal, par_env)
    call set_nnz(mat_sizes, 5)
    call create_matrix(mat_sizes, M)
    call create_vector(vec_sizes, b)
    call create_matrix(mat_sizes, M_exact)
    call create_vector(vec_sizes, b_exact)

    call compute_fluxes(u, v, cell_mesh, BCs, cps, M, b)

    call update(M)
    call update(b)

    call compute_exact_matrix(cell_mesh, flow_direction, discretisation, M_exact, b_exact)

    call update(M_exact)
    call update(b_exact)

    call axpy(-1.0_accs_real, M_exact, M)
    error = norm(M, 1)

    if (error .ge. 1.0e-16) then
      write(message, *) 'FAIL: matrix difference norm too large ', error
      call stop_test(message)
    end if
    
    call axpy(-1.0_accs_real, b_exact, b)
    error = norm(b, 2)

    if (error .ge. 1.0e-16) then
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
  !> @param[out] M             - The resulting matrix
  !> @param[out] b             - The resulting RHS vector
  subroutine compute_exact_matrix(cell_mesh, flow, discretisation, M, b)
    use constants, only: add_mode
    class(mesh), intent(in) :: cell_mesh
    integer(accs_int), intent(in) :: flow
    integer(accs_int), intent(in) :: discretisation
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
    call set_global_size(mat_sizes, cell_mesh%nglobal, cell_mesh%nglobal, par_env)
    call set_nnz(mat_sizes, 5)
    call create_matrix(mat_sizes, M)
    
    call initialise(vec_sizes)
    call set_global_size(vec_sizes, cell_mesh%nglobal, par_env)
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
    if (flow == 1 .and. discretisation == 1) then
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
    else if (flow == 2 .and. discretisation == 1) then
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
    else if (flow == 1 .and. discretisation == 2) then
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
    else if (flow == 2 .and. discretisation == 2) then
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

    diff_coeff = -0.02_accs_real
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
    if (discretisation == 1) then
      adv_coeff = -0.2_accs_real
    else
      adv_coeff = 0.0_accs_real
    endif 

    if (flow == 1) then
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
    diff_coeff = 0.02_accs_real
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
