!> Test that the flux matrix has been computed correctly
!
!> Compares the matrix and RHS calculated for a specified mass flux using the central differencing scheme to the expected solution
program test_compute_fluxes
#include "ccs_macros.inc"

  use testing_lib
  use types, only: field, central_field, face_field
  use mesh_utils, only : build_square_mesh
  use fv, only: compute_fluxes
  use utils, only : update, initialise, finalise, &
       set_size, pack_entries, set_values, zero, &
       set_mode, set_entry, set_row, clear_entries
  use vec, only : create_vector, get_vector_data, restore_vector_data, set_vector_location, create_vector_values
  use mat, only : create_matrix, set_nnz
  use solver, only : axpy, norm
  use constants, only: add_mode, face
  use bc_constants
  use meshing, only: get_global_index, get_local_index, get_boundary_status, &
                     set_cell_location, set_neighbour_location, count_neighbours

  implicit none

  real(ccs_real), parameter :: diffusion_factor = 1.e-2_ccs_real ! XXX: temporarily hard-coded
  type(ccs_mesh) :: mesh
  type(bc_config) :: bcs
  type(vector_spec) :: vec_properties
  class(field), allocatable :: scalar
  class(field), allocatable :: mf
  integer(ccs_int), parameter :: cps = 5
  integer(ccs_int) :: i
  real(ccs_real), dimension(3) :: mf_values

  call init()

  mesh = build_square_mesh(par_env, cps, 1.0_ccs_real)

  bcs%region(1) = bc_region_left
  bcs%region(2) = bc_region_right
  bcs%region(3) = bc_region_top
  bcs%region(4) = bc_region_bottom
  bcs%bc_type(:) = 0
  bcs%bc_type(3) = 1
  bcs%endpoints(:,:) = 1.0_ccs_real
    
  allocate(central_field :: scalar)
  allocate(face_field :: mf)

  call initialise(vec_properties)
  call set_size(par_env, mesh, vec_properties)
  call create_vector(vec_properties, scalar%values)

  call set_vector_location(face, vec_properties)
  call set_size(par_env, mesh, vec_properties)
  call create_vector(vec_properties, mf%values)
  call update(mf%values)

  mf_values = (/ -1.0_ccs_real, 0.0_ccs_real, 1.0_ccs_real /)

  do i = 1, size(mf_values)
    call set_mass_flux(mf, mf_values(i))
    call run_compute_fluxes_test(scalar, mf, mf_values(i), bcs, mesh, cps)
  enddo

  deallocate(scalar)
  deallocate(mf)

  call fin()

  contains

  !> Sets the mass flux array
  subroutine set_mass_flux(mf, mf_value)
    class(field), intent(inout) :: mf     !< The mass flux  
    real(ccs_real) :: mf_value            !< The value to set the mass flux field to
    real(ccs_real), dimension(:), pointer :: mf_data

    call get_vector_data(mf%values, mf_data)
    mf_data(:) = mf_value
    call restore_vector_data(mf%values, mf_data)
    call update(mf%values)
  end subroutine set_mass_flux

  !> Compares the matrix computed for a given velocity field and discretisation to the known solution
  subroutine run_compute_fluxes_test(scalar, mf, mf_value, bcs, mesh, cps)
    class(field), intent(in) :: scalar      !< The scalar field structure
    class(field), intent(in) :: mf          !< The mass flux field
    real(ccs_real), intent(in) :: mf_value  !< The constant value of the mass flux
    class(bc_config), intent(in) :: bcs     !< The BC structure
    type(ccs_mesh), intent(in) :: mesh      !< The mesh structure
    integer(ccs_int), intent(in) :: cps     !< The number of cells per side in the (square) mesh 

    class(ccs_matrix), allocatable :: M, M_exact
    class(ccs_vector), allocatable :: b, b_exact
    type(vector_spec) :: vec_properties
    type(matrix_spec) :: mat_properties
    real(ccs_real) :: error

    
    call initialise(mat_properties)
    call initialise(vec_properties)
    call set_size(par_env, mesh, mat_properties)
    call set_size(par_env, mesh, vec_properties)
    call set_nnz(5, mat_properties)
    call create_matrix(mat_properties, M)
    call create_vector(vec_properties, b)
    call create_matrix(mat_properties, M_exact)
    call create_vector(vec_properties, b_exact)

    call zero(M)

    call compute_fluxes(scalar, mf, mesh, bcs, cps, M, b)

    call update(M)
    call update(b)

    call finalise(M)

    call compute_exact_matrix(mesh, mf_value, cps, M_exact)
    call compute_exact_vector(mesh, mf_value, cps, bcs, b_exact)

    call update(M_exact)
    call update(b_exact)

    call finalise(M_exact)

    call axpy(-1.0_ccs_real, M_exact, M)
    error = norm(M, 1)

    if (error .ge. eps) then
      write(message, *) 'FAIL: matrix difference norm too large ', error
      call stop_test(message)
    end if
    
    call axpy(-1.0_ccs_real, b_exact, b)
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

  !> Computes the expected matrix for a mass flux of 1
  subroutine compute_exact_matrix(mesh, mf_value, cps, M)
    type(ccs_mesh), intent(in) :: mesh      !< The mesh structure
    real(ccs_real), intent(in) :: mf_value  !< The value of the mass flux field
    integer(ccs_int), intent(in) :: cps     !< The number of cells per side
    class(ccs_matrix), intent(inout) :: M   !< The matrix

    ! Local variables
    type(matrix_values) :: mat_values
    type(cell_locator) :: loc_p
    type(neighbour_locator) loc_nb
    integer(ccs_int) :: index_p, index_nb, j, nnb
    integer(ccs_int) :: global_index_p, global_index_nb
    integer(ccs_int) :: sgn
    real(ccs_real) :: face_area, dx
    real(ccs_real) :: adv_coeff, diff_coeff
    real(ccs_real) :: adv_coeff_total, diff_coeff_total
    logical :: is_boundary

    allocate(mat_values%global_row_indices(1))
    allocate(mat_values%global_col_indices(1))
    allocate(mat_values%values(1))

    mat_values%setter_mode = add_mode

    face_area = 1.0_ccs_real/cps

    do index_p = 1, mesh%nlocal
      adv_coeff_total = 0.0_ccs_real
      diff_coeff_total = 0.0_ccs_real
      call set_cell_location(mesh, index_p, loc_p)
      call get_global_index(loc_p, global_index_p)
      call count_neighbours(loc_p, nnb)
      do j = 1, nnb
        call set_neighbour_location(loc_p, j, loc_nb)
        call get_boundary_status(loc_nb, is_boundary)
        if (.not. is_boundary) then
          dx = 1.0_ccs_real/cps
          diff_coeff = -face_area * diffusion_factor / dx
          call get_global_index(loc_nb, global_index_nb)
          call get_local_index(loc_nb, index_nb)
          if (index_nb < index_p) then
            sgn = -1
          else
            sgn = 1
          endif
          adv_coeff = 0.5_ccs_real*mf_value*sgn*face_area
          call pack_entries(1, 1, global_index_p, global_index_nb, adv_coeff + diff_coeff, mat_values)
          call set_values(mat_values, M)
          adv_coeff_total = adv_coeff_total + adv_coeff
          diff_coeff_total = diff_coeff_total + diff_coeff
        else
          dx = 1.0_ccs_real/cps
          diff_coeff = -face_area * diffusion_factor / (0.5_ccs_real * dx)
          adv_coeff = mf_value*face_area
          call pack_entries(1, 1, global_index_p, global_index_p, -(adv_coeff + diff_coeff), mat_values)
          call set_values(mat_values, M)
        endif
      end do
      call pack_entries(1, 1, global_index_p, global_index_p, -(adv_coeff_total + diff_coeff_total), mat_values)
      call set_values(mat_values, M)
    end do
  end subroutine compute_exact_matrix
  
  !> Computes the expected RHS for a mass flux of 1
  subroutine compute_exact_vector(mesh, mf_value, cps, bcs, b)
    type(ccs_mesh), intent(in) :: mesh      !< The mesh structure
    real(ccs_real), intent(in) :: mf_value  !< The value of the mass flux field
    integer(ccs_int), intent(in) :: cps     !< The number of cells per side
    type(bc_config), intent(in) :: bcs      !< The bc structure
    class(ccs_vector), intent(inout) :: b   !< The RHS vector

    type(vector_values) :: vec_values
    type(cell_locator) :: loc_p
    type(neighbour_locator) loc_nb
    integer(ccs_int) :: index_p, j, nnb
    integer(ccs_int) :: global_index_p
    real(ccs_real) :: face_area
    real(ccs_real) :: dx
    real(ccs_real) :: bc_value
    real(ccs_real) :: adv_coeff, diff_coeff
    real(ccs_real) :: adv_coeff_total, diff_coeff_total
    logical :: is_boundary

    call create_vector_values(1_ccs_int, vec_values)
    call set_mode(add_mode, vec_values)

    face_area = 1.0_ccs_real/cps
    do index_p = 1, mesh%nlocal
      call clear_entries(vec_values)
      
      adv_coeff_total = 0.0_ccs_real
      diff_coeff_total = 0.0_ccs_real
      call set_cell_location(mesh, index_p, loc_p)
      call get_global_index(loc_p, global_index_p)
      call count_neighbours(loc_p, nnb)
      do j = 1, nnb
        call set_neighbour_location(loc_p, j, loc_nb)
        call get_boundary_status(loc_nb, is_boundary)
        if (is_boundary) then
          dx = 1.0_ccs_real/cps
          diff_coeff = -face_area * diffusion_factor / (0.5_ccs_real * dx)
          adv_coeff = mf_value*face_area
          if (bcs%bc_type(j) == 0) then
            bc_value = 0.0_ccs_real
          else
            bc_value = 1.0_ccs_real
          endif

          call set_row(global_index_p, vec_values)
          call set_entry(-(adv_coeff + diff_coeff) * bc_value, vec_values)
          call set_values(vec_values, b)
        endif 
      end do
    end do
  end subroutine compute_exact_vector
end program test_compute_fluxes
