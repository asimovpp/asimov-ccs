!> @brief Test that the flux matrix has been computed correctly
!
!> @description Compares the matrix calculated for flows in the +x and +y directions with
!> central and upwind differencing to the known matrix
program test_compute_fluxes
#include "ccs_macros.inc"

  use testing_lib
  use types, only: field, central_field, face_field
  use mesh_utils, only : build_square_mesh
  use fv, only: compute_fluxes, calc_cell_coords
  use utils, only : update, initialise, finalise, &
                set_size, pack_entries, set_values, debug_print, str, zero
  use vec, only : create_vector, get_vector_data, restore_vector_data, set_vector_location
  use mat, only : create_matrix, set_nnz, mat_view
  use solver, only : axpy, norm
  use constants, only: add_mode, insert_mode, face
  use bc_constants
  use petsctypes, only: matrix_petsc
  use meshing, only: get_global_index, get_local_index, get_boundary_status, &
                     set_cell_location, set_neighbour_location, count_neighbours

  implicit none

  real(ccs_real), parameter :: diffusion_factor = 1.e-2_ccs_real ! XXX: temporarily hard-coded
  type(ccs_mesh) :: mesh
  type(bc_config) :: bcs
  type(vector_spec) :: vec_properties
  class(field), allocatable :: scalar
  class(field), allocatable :: u, v
  class(field), allocatable :: mf
  integer(ccs_int), parameter :: cps = 5
  integer(ccs_int) :: direction, discretisation
  integer, parameter :: x_dir = 1, y_dir = 2
  integer, parameter :: central = -1

  call init()

  mesh = build_square_mesh(par_env, cps, 1.0_ccs_real)

  bcs%region(1) = bc_region_left
  bcs%region(2) = bc_region_right
  bcs%region(3) = bc_region_top
  bcs%region(4) = bc_region_bottom
  bcs%bc_type(:) = 0
  bcs%bc_type(3) = 1
  bcs%endpoints(:,:) = 1.0_ccs_real
    
  do direction = x_dir, y_dir
    call dprint("flow direction " // str(direction, "(I3)"))
    discretisation = central
      
    if (discretisation == central) then
      allocate(central_field :: scalar)
      allocate(central_field :: u)
      allocate(central_field :: v)
      allocate(face_field :: mf)
    else
      write(message, *) 'Invalid discretisation type selected'
      call stop_test(message)
    end if

    call initialise(vec_properties)
    call set_size(par_env, mesh, vec_properties)
    call create_vector(vec_properties, scalar%values)
    call create_vector(vec_properties, u%values)
    call create_vector(vec_properties, v%values)

    call set_vector_location(face, vec_properties)
    call set_size(par_env, mesh, vec_properties)
    call create_vector(vec_properties, mf%values)
    call update(mf%values)

    call set_fields(mesh, direction, u, v, mf)
    call run_compute_fluxes_test(scalar, u, v, mf, bcs, mesh, cps, direction, discretisation)

    deallocate(scalar)
    deallocate(u)
    deallocate(v)
    deallocate(mf)
  end do

  call fin()

  contains

  !> @brief Sets the velocity field in the desired direction and discretisation
  !
  !> @param[in] mesh - The mesh structure
  !> @param[in] direction - Integer indicating the direction of the velocity field
  !> @param[out] u, v     - The velocity fields in x and y directions
  !> @param[out] mf       - The mass flux field 
  subroutine set_fields(mesh, direction, u, v, mf)
    use meshing, only: set_cell_location, get_global_index
    class(ccs_mesh), intent(in) :: mesh
    integer(ccs_int), intent(in) :: direction
    class(field), intent(inout) :: u, v, mf
    type(cell_locator) :: loc_p
    type(vector_values) :: u_vals, v_vals, mf_vals
    integer(ccs_int) :: index_p, global_index_p
    real(ccs_real) :: u_val, v_val, mf_val
    real(ccs_real), dimension(:), pointer :: mf_data

    u_vals%setter_mode = insert_mode
    v_vals%setter_mode = insert_mode
    
    associate(n_local => mesh%nlocal)
      allocate(u_vals%indices(n_local))
      allocate(v_vals%indices(n_local))
      allocate(u_vals%values(n_local))
      allocate(v_vals%values(n_local))
      
      ! Set IC velocity fields
      do index_p = 1, n_local
        call set_cell_location(mesh, index_p, loc_p)
        call get_global_index(loc_p, global_index_p)

        if (direction == x_dir) then
          u_val = 1.0_ccs_real
          v_val = 0.0_ccs_real
        else if (direction == y_dir) then
          u_val = 0.0_ccs_real
          v_val = 1.0_ccs_real
        end if

        !u_val = 0.0_ccs_real
        !v_val = 0.0_ccs_real
        
        call pack_entries(index_p, global_index_p, u_val, u_vals)
        call pack_entries(index_p, global_index_p, v_val, v_vals)
      end do
    end associate
    call set_values(u_vals, u%values)
    call set_values(v_vals, v%values)

    call get_vector_data(mf%values, mf_data)
    mf_data(:) = 1.0_ccs_real
    call restore_vector_data(mf%values, mf_data)

    call update(u%values)
    call update(v%values)
    call update(mf%values)
    
    deallocate(u_vals%indices, v_vals%indices, u_vals%values, v_vals%values)
  end subroutine set_fields

  !> @brief Compares the matrix computed for a given velocity field and discretisation to the known solution
  !
  !> @param[in] scalar         - The scalar field structure
  !> @param[in] u, v           - The velocity field structures
  !> @param[in] bcs            - The BC structure
  !> @param[in] mesh      - The mesh structure
  !> @param[in] cps            - The number of cells per side in the (square) mesh 
  !> @param[in] flow_direction - Integer indicating the direction of the flow 
  !> @param[in] discretisation - Integer indicating the discretisation scheme being tested 
  subroutine run_compute_fluxes_test(scalar, u, v, mf, bcs, mesh, cps, flow_direction, discretisation)
    class(field), intent(in) :: scalar
    class(field), intent(in) :: u, v
    class(field), intent(in) :: mf
    class(bc_config), intent(in) :: bcs
    type(ccs_mesh), intent(in) :: mesh
    integer(ccs_int), intent(in) :: cps
    integer(ccs_int), intent(in) :: flow_direction
    integer(ccs_int), intent(in) :: discretisation

    class(ccs_matrix), allocatable :: M, M_exact
    class(ccs_vector), allocatable :: b, b_exact
    type(vector_spec) :: vec_properties
    type(matrix_spec) :: mat_properties
    real(ccs_real) :: error
    logical :: assembled
    integer(ccs_int), dimension(2) :: rows, cols
    integer(ccs_int), dimension(2) :: sub_M, sub_M_exact

    
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

    call compute_exact_matrix(mesh, flow_direction, discretisation, cps, M_exact)
    call compute_exact_vector(mesh, flow_direction, discretisation, cps, bcs, b_exact)

    call update(M_exact)
    call update(b_exact)

    call finalise(M_exact)

    call mat_view(M_exact, "M_exact.dat")
    call mat_view(M, "M.dat")

    call axpy(-1.0_ccs_real, M_exact, M)
    error = norm(M, 1)
    call dprint("norm " // str(error))

    if (error .ge. eps) then
      write(message, *) 'FAIL: matrix difference norm too large ', error
      call stop_test(message)
    end if
    
    call axpy(-1.0_ccs_real, b_exact, b)
    error = norm(b, 2)
    call dprint("norm " // str(error))

    if (error .ge. eps) then
      write(message, *) 'FAIL: vector difference norm too large ', error
      call stop_test(message)
    end if

    deallocate(M)
    deallocate(b)
    deallocate(M_exact)
    deallocate(b_exact)
  end subroutine run_compute_fluxes_test

  subroutine compute_exact_matrix(mesh, flow, discretisation, cps, M)
    type(ccs_mesh), intent(in) :: mesh
    integer(ccs_int), intent(in) :: flow
    integer(ccs_int), intent(in) :: discretisation
    integer(ccs_int), intent(in) :: cps
    class(ccs_matrix), intent(inout) :: M

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

    allocate(mat_values%row_indices(1))
    allocate(mat_values%col_indices(1))
    allocate(mat_values%values(1))

    mat_values%setter_mode = add_mode

    face_area = 1.0_ccs_real/cps

    do index_p = 1, cps**2
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
          if (index_p > index_nb) then
            sgn = 1
          else
            sgn = -1
          endif
          adv_coeff = -0.5*sgn*face_area
          call pack_entries(1, 1, global_index_p, global_index_nb, adv_coeff + diff_coeff, mat_values)
          call set_values(mat_values, M)
          adv_coeff_total = adv_coeff_total + adv_coeff
          diff_coeff_total = diff_coeff_total + diff_coeff
        else
          dx = 1.0_ccs_real/cps
          diff_coeff = -face_area * diffusion_factor / (0.5_ccs_real * dx)
          adv_coeff = face_area
          call pack_entries(1, 1, global_index_p, global_index_p, -(adv_coeff + diff_coeff), mat_values)
          call set_values(mat_values, M)
        endif
      end do
      call pack_entries(1, 1, global_index_p, global_index_p, -(adv_coeff_total + diff_coeff_total), mat_values)
      call set_values(mat_values, M)
    end do
  end subroutine compute_exact_matrix
  
  subroutine compute_exact_vector(mesh, flow, discretisation, cps, bcs, b)
    type(ccs_mesh), intent(in) :: mesh
    integer(ccs_int), intent(in) :: flow
    integer(ccs_int), intent(in) :: discretisation
    integer(ccs_int), intent(in) :: cps
    type(bc_config), intent(in) :: bcs
    class(ccs_vector), intent(inout) :: b

    type(vector_values) :: vec_values
    type(cell_locator) :: loc_p
    type(neighbour_locator) loc_nb
    integer(ccs_int) :: index_p, index_nb, j, nnb
    integer(ccs_int) :: global_index_p, global_index_nb
    real(ccs_real) :: face_area
    real(ccs_real) :: dx
    real(ccs_real) :: bc_value
    real(ccs_real) :: adv_coeff, diff_coeff
    real(ccs_real) :: adv_coeff_total, diff_coeff_total
    logical :: is_boundary

    allocate(vec_values%indices(1))
    allocate(vec_values%values(1))

    vec_values%setter_mode = add_mode

    face_area = 1.0_ccs_real/cps
    do index_p = 1, cps**2
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
          adv_coeff = face_area
          if (bcs%bc_type(j) == 0) then
            bc_value = 0.0_ccs_real
          else
            bc_value = 1.0_ccs_real
          endif
          call pack_entries(1, global_index_p, -(adv_coeff + diff_coeff)*bc_value, vec_values)
          call set_values(vec_values, b)
        endif 
      end do
    end do
  end subroutine compute_exact_vector
end program test_compute_fluxes
