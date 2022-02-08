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

  local_idx = int(0.5*square_mesh%nlocal + 2, accs_int)
  do direction = 1, 2
    do discretisation = 1, 2
      call set_velocity_fields(direction, cps, discretisation, u, v)
      do ngb = 1, 4
        call set_cell_indices(local_idx, ngb, self_idx, ngb_idx, face_surface_area, normal)
        call calc_cell_coords(self_idx, cps, row, col)
        call run_advection_coeff_test(u, v, self_idx, row, col, ngb_idx, face_surface_area, normal)
      end do
      call tidy_velocity_fields(u, v)
    end do
  end do

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

  subroutine run_advection_coeff_test(u, v, self_idx, row, col, ngb_idx, face_surface_area, normal)
    class(field), intent(in) :: u, v
    integer(accs_int), intent(in) :: self_idx
    integer(accs_int), intent(in) :: row, col
    integer(accs_int), intent(in) :: ngb_idx
    real(accs_real), intent(in) :: face_surface_area
    real(accs_real), dimension(ndim), intent(in) :: normal

    real(accs_real) :: coeff
    real(accs_real) :: expected_coeff

    select type(u)
      type is(central_field)
        select type(v)
          type is(central_field)
            call calc_advection_coeff(self_idx, ngb_idx, face_surface_area, cps, u, v, 0_accs_int, coeff)
          class default
            write(message, *) "FAIL: incorrect velocity field discretisation"
            call stop_test(message)
        end select
      type is(upwind_field)
        select type(v)
          type is(upwind_field)
            call calc_advection_coeff(self_idx, ngb_idx, face_surface_area, cps, u, v, 0_accs_int, coeff)
          class default
            write(message, *) "FAIL: incorrect velocity field discretisation"
            call stop_test(message)
        end select
      class default
        write(message, *) "FAIL: incorrect velocity field discretisation"
        call stop_test(message)
    end select

    select type(u)
      type is(upwind_field)
        expected_coeff = min(-face_surface_area*(u%val(col,row)*normal(1) + v%val(col,row)*normal(2)), 0.0_accs_real)
      type is(central_field)
        expected_coeff = -face_surface_area*0.5_accs_real*(u%val(col,row)*normal(1) + v%val(col,row)*normal(2))
      class default
        write(message, *) "FAIL: incorrect velocity field discretisation"
        call stop_test(message)
    end select
        
    select type(u)
      type is(upwind_field)
        if (abs(coeff - expected_coeff) .ge. tiny(coeff)) then
          write(message, *) "FAIL: incorrect upwind advection coefficient computed. Expected ", expected_coeff, " computed ", &
                            coeff, " normal ", normal, " u ", u%val(col, row), " v ", v%val(col, row)
          call stop_test(message)
        end if
      type is(central_field)
        if (abs(coeff - expected_coeff) .ge. tiny(coeff)) then
          write(message, *) "FAIL: incorrect central advection coefficient computed. Expected ", expected_coeff, " computed ", &
                            coeff, " normal ", normal, " u ", u%val(col, row), " v ", v%val(col, row)
          call stop_test(message)
        end if
      class default
        write(message, *) "FAIL: incorrect velocity field discretisation"
        call stop_test(message)
    end select
        

  end subroutine run_advection_coeff_test

end program test_mesh_neighbours
