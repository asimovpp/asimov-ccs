!> @brief Test that cells have correct numbers of neighbours
!
!> @description for any mesh with >1 cell, every cell must have at least 1 neighbour.
program test_compute_fluxes

  use testing_lib
  use constants, only: ndim
  use types, only: field, upwind_field, central_field, cell_locator, face_locator, neighbour_locator, &
                   set_cell_location, set_face_location, set_neighbour_location
  use mesh_utils, only : build_square_mesh, global_index, face_area, face_normal
  use fv, only: compute_fluxes, calc_cell_coords
  use BC_constants
  use utils, only : update, finalise, initialise, &
                set_global_size, pack_entries, set_values, axpy, norm
  use vec, only : create_vector
  use mat, only : create_matrix, set_nnz


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

  subroutine tidy_velocity_fields(u, v)
    class(field), allocatable :: u, v

    deallocate(u%val)
    deallocate(v%val)
    deallocate(u)
    deallocate(v)
  end subroutine tidy_velocity_fields

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
    integer(accs_int) :: self_idx, ngb_idx, local_idx
    integer(accs_int) :: ngb
    integer(accs_int) :: self_row, self_col
    integer(accs_int) :: ngb_row, ngb_col
    real(accs_real) :: face_surface_area
    real(accs_real) :: mat_val
    real(accs_real) :: mat_expected
    real(accs_real) :: error
    real(accs_real), dimension(ndim) :: normal
    
    call initialise(mat_sizes)
    call initialise(vec_sizes)
    call set_global_size(mat_sizes, cell_mesh%n, cell_mesh%n, par_env)
    call set_global_size(vec_sizes, cell_mesh%n, par_env)
    call set_nnz(mat_sizes, 5)
    call create_matrix(mat_sizes, M)
    call create_vector(vec_sizes, b)
    call create_matrix(mat_sizes, M_exact)
    call create_vector(vec_sizes, b_exact)

    call compute_fluxes(u, v, cell_mesh, BCs, cps, M, b)

    call update(M)
    call update(b)

    call compute_exact_matrix(cell_mesh, flow_direction, discretisation, M_exact)

    call update(M_exact)
    call update(b_exact)

    call axpy(-1.0_accs_real, M_exact, M)
    error = norm(M, 1)

    if (error .ge. 1.0e-16) then
      write(message, *) 'FAIL: matrix difference norm too large ', error
      call stop_test(message)
    end if

    deallocate(M)
    deallocate(b)
    deallocate(M_exact)
    deallocate(b_exact)
  end subroutine run_compute_fluxes_test

  subroutine compute_exact_matrix(cell_mesh, flow, discretisation, M)
    use constants, only: add_mode
    class(mesh), intent(in) :: cell_mesh
    integer(accs_int), intent(in) :: flow
    integer(accs_int), intent(in) :: discretisation
    class(matrix), allocatable :: M

    type(matrix_init_data) :: mat_sizes
    type(matrix_values) :: mat_coeffs
    type(cell_locator) :: self_loc
    type(neighbour_locator) :: ngb_loc
    real(accs_real) :: diff_coeff
    real(accs_real) :: proc_zero_factor ! set to zero for non-zero ranks so that we're not double counting when running in parallel
    integer(accs_int) :: i, j
    integer(accs_int) :: global_idx, ngb_idx
    integer(accs_int) :: mat_counter

    call initialise(mat_sizes)
    call set_global_size(mat_sizes, cell_mesh%n, cell_mesh%n, par_env)
    call set_nnz(mat_sizes, 5)
    call create_matrix(mat_sizes, M)

    mat_coeffs%mode = add_mode
    
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
        do i = 1, cell_mesh%n
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
        do i = 1, cell_mesh%n
          mat_counter = 1
          if (i .le. cps) then
            call pack_entries(mat_coeffs, 1, mat_counter, i, i, -0.3_accs_real*proc_zero_factor)
            mat_counter = mat_counter + 1
          else if (i > cell_mesh%n - cps) then
            call pack_entries(mat_coeffs, 1, mat_counter, i, i, -0.1_accs_real*proc_zero_factor)
            mat_counter = mat_counter + 1
          end if

          if (i + cps .le. cell_mesh%n) then
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
        do i = 1, cell_mesh%n
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
        do i = cps+1, cell_mesh%n
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

    diff_coeff = -0.02
    ! Diffusion coefficients
    do i = 1, cell_mesh%n
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

      if (i + 1 .le. cell_mesh%n .and. mod(i, cps) .ne. 0) then
        call pack_entries(mat_coeffs, 1, mat_counter, i, i+1, diff_coeff*proc_zero_factor)
        mat_counter = mat_counter + 1
      end if
      if (i + cps .le. cell_mesh%n) then
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
    
  end subroutine compute_exact_matrix

end program test_compute_fluxes
