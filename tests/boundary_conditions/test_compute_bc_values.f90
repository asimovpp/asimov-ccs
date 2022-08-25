program test_compute_bc_values
#include "ccs_macros.inc"

  use testing_lib
  use kinds, only: ccs_real, ccs_int
  use types, only : field, central_field, face_locator, cell_locator, neighbour_locator, vector_spec
  use utils, only : debug_print, exit_print, str, update, set_size, initialise
  use fv, only: compute_boundary_values
  use constants, only: ndim, cell
  use bc_constants
  use boundary_conditions, only: allocate_bc_arrays
  use meshing, only: get_local_index, set_cell_location, count_neighbours, set_neighbour_location, get_boundary_status, set_face_location, get_face_normal
  use mesh_utils, only: build_square_mesh
  use vec, only: create_vector, set_vector_location, get_vector_data, restore_vector_data

  class(field), allocatable :: dirichlet_field, extrapolated_field, sym_field
  integer(ccs_int) :: component
  integer(ccs_int) :: index_p, index_nb
  integer(ccs_int) :: nnb
  integer(ccs_int), parameter :: cps = 10, n_boundaries = 4
  integer(ccs_int) :: j
  type(face_locator) :: loc_f
  type(cell_locator) :: loc_p
  type(neighbour_locator) :: loc_nb
  type(ccs_mesh) :: mesh
  type(vector_spec) :: vec_properties
  real(ccs_real), dimension(ndim) :: face_norm
  real(ccs_real) :: bc_val, expected_bc_value
  real(ccs_real), dimension(:), pointer :: extrapolated_field_data
  real(ccs_real), dimension(:), pointer :: x_gradient_data
  real(ccs_real), dimension(:), pointer :: y_gradient_data
  real(ccs_real), dimension(:), pointer :: z_gradient_data
  real(ccs_real), dimension(:), pointer :: sym_field_data
  logical :: is_boundary

  call init()

  mesh = build_square_mesh(par_env, cps, 1.0_ccs_real)

  ! set locations
  index_p = 3
  call set_cell_location(mesh, index_p, loc_p)
  call count_neighbours(loc_p, nnb)
  do j = 1, nnb
    call set_neighbour_location(loc_p, j, loc_nb)
    call get_boundary_status(loc_nb, is_boundary)
    if (is_boundary) then
      call set_face_location(mesh, index_p, j, loc_f)
      call get_face_normal(loc_f, face_norm)
      call get_local_index(loc_nb, index_nb)
	    exit
	  end if
  end do 

  ! Check Dirichlet BC
  component = 1
  expected_bc_value = 7.5
  allocate(central_field :: dirichlet_field)
  call allocate_bc_arrays(n_boundaries, dirichlet_field%bcs)
  dirichlet_field%bcs%bc_types = bc_type_dirichlet
  dirichlet_field%bcs%values = expected_bc_value
  dirichlet_field%bcs%ids = (/ (j, j = 1, n_boundaries) /)

  call compute_boundary_values(dirichlet_field, component, index_nb, index_p, loc_p, loc_f, face_norm, bc_val)
  call assert_equal(bc_val, expected_bc_value, '("bc values do not match received ", f7.4, " expected ", f7.4)')
  call dprint("done dirichlet test")

  ! Check extrapolated BC
  call initialise(vec_properties)
  call set_vector_location(cell, vec_properties)
  call set_size(par_env, mesh, vec_properties)
  allocate(central_field :: extrapolated_field)
  call create_vector(vec_properties, extrapolated_field%values)
  call create_vector(vec_properties, extrapolated_field%x_gradients)
  call create_vector(vec_properties, extrapolated_field%y_gradients)
  call create_vector(vec_properties, extrapolated_field%z_gradients)
  call update(extrapolated_field%values)
  call update(extrapolated_field%x_gradients)
  call update(extrapolated_field%y_gradients)
  call update(extrapolated_field%z_gradients)
  call allocate_bc_arrays(n_boundaries, extrapolated_field%bcs)
  extrapolated_field%bcs%bc_types = bc_type_extrapolate
  extrapolated_field%bcs%ids = (/ (j, j = 1, n_boundaries) /)

  call get_vector_data(extrapolated_field%values, extrapolated_field_data)
  call get_vector_data(extrapolated_field%x_gradients, x_gradient_data)
  call get_vector_data(extrapolated_field%y_gradients, y_gradient_data)
  call get_vector_data(extrapolated_field%z_gradients, z_gradient_data)
  x_gradient_data = 0
  y_gradient_data = 1.0_ccs_real / cps
  z_gradient_data = 0
  do j = 1, mesh%nlocal
    extrapolated_field_data(j) = j / cps
  end do
  expected_bc_value = extrapolated_field_data(index_p) - 0.5_ccs_real / cps * y_gradient_data(index_p)
  call restore_vector_data(extrapolated_field%values, extrapolated_field_data)
  
  call compute_boundary_values(extrapolated_field, component, index_nb, index_p, loc_p, loc_f, face_norm, bc_val, &
                               x_gradient_data, y_gradient_data, z_gradient_data)
  
  call restore_vector_data(extrapolated_field%x_gradients, x_gradient_data)
  call restore_vector_data(extrapolated_field%y_gradients, y_gradient_data)
  call restore_vector_data(extrapolated_field%z_gradients, z_gradient_data)

  call assert_equal(bc_val, expected_bc_value, '("bc values do not match received ", f7.4, " expected ", f7.4)')
  call dprint("done extrapolated test")

  ! Check symmetric test
  allocate(central_field :: sym_field)
  call create_vector(vec_properties, sym_field%values)
  call update(sym_field%values)
  call allocate_bc_arrays(n_boundaries, sym_field%bcs)
  sym_field%bcs%bc_types = bc_type_sym
  sym_field%bcs%ids = (/ (j, j = 1, n_boundaries) /)

  call get_vector_data(sym_field%values, sym_field_data)
  do j = 1, mesh%nlocal
    sym_field_data(j) = j / cps + 1
  end do
  call restore_vector_data(sym_field%values, sym_field_data)
  
  do j = 1, ndim
    component = j
    if (j == 2) then
      expected_bc_value = 0
    else 
      expected_bc_value = 1
    end if
    call compute_boundary_values(sym_field, component, index_nb, index_p, loc_p, loc_f, face_norm, bc_val)
    call assert_equal(bc_val, expected_bc_value, '("bc values do not match received ", f7.4, " expected ", f7.4)')
  end do
  call dprint("done symmetric test")


  call dprint("done")

  call fin()

end program test_compute_bc_values