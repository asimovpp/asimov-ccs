!> @brief Test that cells have correct numbers of neighbours
!
!> @description for any mesh with >1 cell, every cell must have at least 1 neighbour.
program test_mesh_neighbours

  use testing_lib
  use constants, only: ndim
  use types, only: field, upwind_field, central_field, cell_locator, face_locator, neighbour_locator, &
                   set_cell_location, set_face_location, set_neighbour_location
  use mesh_utils, only : build_square_mesh, global_index, face_area, face_normal
  use fv, only: calc_advection_coeff, calc_cell_coords

  type(mesh) :: square_mesh
  class(field), allocatable :: u, v
  integer(accs_int), parameter :: cps = 50
  integer(accs_int) :: self_idx, ngb_idx, local_idx
  integer(accs_int) :: ngb
  integer(accs_int) :: direction, discretisation
  integer(accs_int) :: row, col
  real(accs_real) :: face_surface_area
  real(accs_real), dimension(ndim) :: normal

  call init()

  square_mesh = build_square_mesh(cps, 1.0_accs_real, par_env)


  call fin()

  contains

  subroutine set_cell_indices(local_idx, ngb, self_idx, ngb_idx, face_surface_area, normal)
    integer(accs_int), intent(in) :: local_idx
    integer(accs_int), intent(in) :: ngb
    integer(accs_int), intent(out) :: self_idx
    integer(accs_int), intent(out) :: ngb_idx
    real(accs_real), intent(out) :: face_surface_area
    real(accs_real), intent(out), dimension(ndim) :: normal

    type(cell_locator) :: self_loc
    type(neighbour_locator) :: ngb_loc
    type(face_locator) :: face_loc
    
    call set_cell_location(self_loc, square_mesh, local_idx)
    call global_index(self_loc, self_idx)
  
    call set_neighbour_location(ngb_loc, self_loc, ngb)
    call global_index(ngb_loc, ngb_idx)

    call set_face_location(face_loc, square_mesh, local_idx, ngb)
    call face_area(face_loc, face_surface_area)

    call face_normal(face_loc, normal)
  end subroutine set_cell_indices

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

  subroutine run_compute_fluxes_test(u, v, cell_mesh, cps, M, b)
    class(field), intent(in) :: u, v
    type(mesh), intent(in) :: cell_mesh
    integer(accs_int), intent(in) :: cps
    class(matrix), intent(in) :: M
    class(vector), intent(in) :: b

    call compute_fluxes(u, v, cell_mesh, cps, M, b)

    call update(M)
    call update(b)

    call set_exact_solution(M_exact)
    call set_exact_solution(b_exact)

    mat_residual = norm(M - M_exact, 1)
    vec_residual = norm(b - b_exact, 2)
    if (mat_residual .ge. tiny(1.0_accs_real)) then
      write(message, *) "FAIL: incorrect matrix computed. residual ", mat_residual
      call stop_test(message)
    end if
    if (vec_residual .ge. tiny(1.0_accs_real)) then
      write(message, *) "FAIL: incorrect vector computed. residual ", vec_residual
      call stop_test(message)
    end if
  end subroutine run_advection_coeff_test

end program test_mesh_neighbours
