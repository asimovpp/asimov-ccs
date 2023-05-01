program test_compute_bc_values
#include "ccs_macros.inc"

  use testing_lib
  use kinds, only: ccs_real, ccs_int
  use types, only: field, central_field, face_locator, cell_locator, neighbour_locator, vector_spec
  use utils, only: debug_print, exit_print, str, update, set_size, initialise
  use fv, only: compute_boundary_values
  use constants, only: ndim, cell
  use bc_constants
  use boundary_conditions, only: allocate_bc_arrays
  use meshing, only: get_local_index, get_cell_location, count_neighbours, set_neighbour_location, &
                     get_boundary_status, get_face_location, get_face_normal, get_local_num_cells
  use mesh_utils, only: build_square_mesh
  use vec, only: create_vector, set_vector_location, get_vector_data, restore_vector_data

  implicit none

  integer(ccs_int) :: local_num_cells
  integer(ccs_int) :: nnb
  integer(ccs_int), parameter :: cps = 10, n_boundaries = 4
  integer(ccs_int) :: j, k
  type(face_locator) :: loc_f
  type(cell_locator) :: loc_p
  type(neighbour_locator) :: loc_nb
  type(ccs_mesh) :: mesh
  real(ccs_real), dimension(ndim) :: face_normal
  logical :: is_boundary

  integer(ccs_int) :: index_p

  call init()

  mesh = build_square_mesh(par_env, cps, 1.0_ccs_real)

  ! set locations
  index_p = 0
  call get_local_num_cells(mesh, local_num_cells)
  do k = 1, local_num_cells
    call get_cell_location(mesh, k, loc_p)
    call count_neighbours(loc_p, nnb)
    do j = 1, nnb
      call set_neighbour_location(loc_p, j, loc_nb)
      call get_boundary_status(loc_nb, is_boundary)
      call get_face_location(mesh, k, j, loc_f)
      call get_face_normal(loc_f, face_normal)
      if (is_boundary .and. all(face_normal .eq. (/0, -1, 0/))) then
        index_p = k
        exit
      end if
    end do
    if (index_p /= 0) then
      exit
    end if
  end do

  if (index_p /= 0) then
    call check_dirichlet_bc(loc_p, loc_f)
    call check_neumann_bc(loc_p, loc_f)
    call check_extrapolated_bc(loc_p, loc_f, cps)
    !call check_symmetric_bc(loc_p, loc_f, cps)  ! XXX: fix symmetric BC implementation and test accordingly

    call dprint("done")
  end if

  call fin()

contains

  ! Checks whether the dirichlet bcs are being computed correctly
  subroutine check_dirichlet_bc(loc_p, loc_f)
    type(cell_locator), intent(in) :: loc_p   !< cell location to check bc at
    type(face_locator), intent(in) :: loc_f   !< the face location at the boundary

    type(vector_spec) :: vec_properties
    integer(ccs_int) :: component = 1
    real(ccs_real) :: expected_bc_value = 7.5
    real(ccs_real) :: bc_val
    class(field), allocatable :: dirichlet_field
    integer(ccs_int), parameter :: n_boundaries = 4
    real(ccs_real), dimension(ndim) :: face_norm

    call get_face_normal(loc_f, face_norm)
    call initialise(vec_properties)
    call set_vector_location(cell, vec_properties)
    call set_size(par_env, mesh, vec_properties)
    allocate (central_field :: dirichlet_field)
    call create_vector(vec_properties, dirichlet_field%values)
    call update(dirichlet_field%values)
    call allocate_bc_arrays(n_boundaries, dirichlet_field%bcs)
    call create_vector(vec_properties, dirichlet_field%values)
    dirichlet_field%bcs%bc_types = bc_type_dirichlet
    dirichlet_field%bcs%values = expected_bc_value
    dirichlet_field%bcs%ids = (/(j, j=1, n_boundaries)/)

    call compute_boundary_values(dirichlet_field, component, loc_p, loc_f, face_norm, bc_val)
    call assert_eq(bc_val, expected_bc_value, "bc values do not match received")
    call dprint("done dirichlet test")
  end subroutine check_dirichlet_bc

  ! Checks whether the neumann bcs are being computed correctly
  subroutine check_neumann_bc(loc_p, loc_f)
    type(cell_locator), intent(in) :: loc_p   !< cell location to check bc at
    type(face_locator), intent(in) :: loc_f   !< the face location at the boundary

    type(vector_spec) :: vec_properties
    integer(ccs_int) :: component = 1
    real(ccs_real) :: expected_bc_value = 7.5
    real(ccs_real) :: bc_val
    class(field), allocatable :: neumann_field
    real(ccs_real), dimension(:), pointer :: neumann_field_data
    integer(ccs_int), parameter :: n_boundaries = 4
    real(ccs_real), dimension(ndim) :: face_norm
    integer(ccs_int) :: j

    call initialise(vec_properties)
    call set_vector_location(cell, vec_properties)
    call set_size(par_env, mesh, vec_properties)
    allocate (central_field :: neumann_field)
    call create_vector(vec_properties, neumann_field%values)
    call update(neumann_field%values)
    call allocate_bc_arrays(n_boundaries, neumann_field%bcs)
    neumann_field%bcs%bc_types = bc_type_neumann
    neumann_field%bcs%values = 0.0_ccs_real
    neumann_field%bcs%ids = (/(j, j=1, n_boundaries)/)

    call get_vector_data(neumann_field%values, neumann_field_data)
    do j = 1, local_num_cells
      neumann_field_data(j) = expected_bc_value
    end do
    call restore_vector_data(neumann_field%values, neumann_field_data)
    call update(neumann_field%values)
    call dprint("set neumann field")

    call compute_boundary_values(neumann_field, component, loc_p, loc_f, face_norm, bc_val)
    call assert_eq(bc_val, expected_bc_value, "bc values do not match received")
    call dprint("done neumann test")
  end subroutine check_neumann_bc

  ! Checks whether extrapolation bcs are being computed correctly
  subroutine check_extrapolated_bc(loc_p, loc_f, cps)
    type(cell_locator), intent(in) :: loc_p   !< cell location to check bc at
    type(face_locator), intent(in) :: loc_f   !< the face location at the boundary
    integer(ccs_int), intent(in) :: cps       !< the number of cells per side of the mesh

    class(field), allocatable :: extrapolated_field
    integer(ccs_int) :: component = 1
    type(vector_spec) :: vec_properties
    real(ccs_real) :: bc_val, expected_bc_value
    real(ccs_real), dimension(:), pointer :: extrapolated_field_data
    real(ccs_real), dimension(:), pointer :: x_gradient_data
    real(ccs_real), dimension(:), pointer :: y_gradient_data
    real(ccs_real), dimension(:), pointer :: z_gradient_data
    integer(ccs_int), parameter :: n_boundaries = 4
    real(ccs_real), dimension(ndim) :: face_norm

    associate (mesh => loc_f%mesh)
      call initialise(vec_properties)
      call set_vector_location(cell, vec_properties)
      call set_size(par_env, mesh, vec_properties)
      allocate (central_field :: extrapolated_field)
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
      extrapolated_field%bcs%ids = (/(j, j=1, n_boundaries)/)

      call get_face_normal(loc_f, face_norm)

      call get_vector_data(extrapolated_field%values, extrapolated_field_data)
      call get_vector_data(extrapolated_field%x_gradients, x_gradient_data)
      call get_vector_data(extrapolated_field%y_gradients, y_gradient_data)
      call get_vector_data(extrapolated_field%z_gradients, z_gradient_data)
      x_gradient_data = 0
      y_gradient_data = 1.0_ccs_real / cps
      z_gradient_data = 0
      do j = 1, local_num_cells
        extrapolated_field_data(j) = j / cps
      end do
      expected_bc_value = extrapolated_field_data(index_p) - 0.5_ccs_real / cps * y_gradient_data(index_p)
      call restore_vector_data(extrapolated_field%values, extrapolated_field_data)

      call compute_boundary_values(extrapolated_field, component, loc_p, loc_f, face_norm, bc_val, &
                                   x_gradient_data, y_gradient_data, z_gradient_data)

      call restore_vector_data(extrapolated_field%x_gradients, x_gradient_data)
      call restore_vector_data(extrapolated_field%y_gradients, y_gradient_data)
      call restore_vector_data(extrapolated_field%z_gradients, z_gradient_data)

      call assert_eq(bc_val, expected_bc_value, "bc values do not match received")
      call dprint("done extrapolated test")
    end associate
  end subroutine check_extrapolated_bc

  ! ! Checks whether symmetric bcs are being computed correctly
  ! subroutine check_symmetric_bc(loc_p, loc_f, cps)
  !   type(cell_locator), intent(in) :: loc_p   !< cell location to check bc at
  !   type(face_locator), intent(in) :: loc_f   !< the face location at the boundary
  !   integer(ccs_int), intent(in) :: cps       !< the number of cells per side of the mesh

  !   type(vector_spec) :: vec_properties
  !   class(field), allocatable :: sym_field
  !   integer(ccs_int) :: component
  !   integer(ccs_int) :: n_boundaries = 4
  !   real(ccs_real), dimension(:), pointer :: sym_field_data
  !   real(ccs_real), dimension(ndim) :: face_norm
  !   real(ccs_real) :: expected_bc_value, bc_val

  !   call get_face_normal(loc_f, face_norm)

  !   associate (mesh => loc_f%mesh)
  !     call initialise(vec_properties)
  !     call set_vector_location(cell, vec_properties)
  !     call set_size(par_env, mesh, vec_properties)
  !     allocate (central_field :: sym_field)
  !     call create_vector(vec_properties, sym_field%values)
  !     call update(sym_field%values)
  !     call allocate_bc_arrays(n_boundaries, sym_field%bcs)
  !     sym_field%bcs%bc_types = bc_type_sym
  !     sym_field%bcs%ids = (/(j, j=1, n_boundaries)/)

  !     call get_vector_data(sym_field%values, sym_field_data)
  !     do j = 1, local_num_cells
  !       sym_field_data(j) = j / cps + 1
  !     end do
  !     call restore_vector_data(sym_field%values, sym_field_data)

  !     do j = 1, ndim
  !       component = j
  !       if (j == 2) then
  !         expected_bc_value = 0
  !       else
  !         expected_bc_value = 1
  !       end if
  !       call compute_boundary_values(sym_field, component, loc_p, loc_f, face_norm, bc_val)
  !       call assert_eq(bc_val, expected_bc_value, "bc values do not match received")
  !     end do
  !     call dprint("done symmetric test")
  !   end associate

  ! end subroutine check_symmetric_bc

end program test_compute_bc_values
