!> @brief Test that advection coefficients are calculated correctly
!
!> @description Computes the advection coefficients for two flow directions
!> (in +x, +y directions) for central and upwind differencing and compares to
!> known values
program test_advection_coeff

  use testing_lib
  use constants, only: ndim
  use types, only: field, upwind_field, central_field, cell_locator, face_locator, neighbour_locator
  use mesh_utils, only : build_square_mesh
  use fv, only: calc_advection_coeff, calc_cell_coords
  use meshing, only: set_cell_location, set_face_location, set_neighbour_location, &
                     get_global_index, get_face_area, get_face_normal

  type(mesh) :: square_mesh
  class(field), allocatable :: u, v
  integer(accs_int), parameter :: cps = 50
  integer(accs_int) :: self_idx, ngb_idx, local_idx
  integer(accs_int) :: ngb
  integer(accs_int) :: direction, discretisation
  integer(accs_int) :: row, col
  real(accs_real) :: face_area
  real(accs_real), dimension(ndim) :: normal
  
  call init()

  square_mesh = build_square_mesh(cps, 1.0_accs_real, par_env)

  local_idx = int(0.5*square_mesh%nlocal + 2, accs_int)
  do direction = 1, 2
    do discretisation = 1, 2
      call set_velocity_fields(direction, cps, discretisation, u, v)
      do ngb = 1, 4
        call set_cell_indices(local_idx, ngb, self_idx, ngb_idx, face_area, normal)
        call calc_cell_coords(self_idx, cps, row, col)
        call run_advection_coeff_test(u, v, self_idx, row, col, ngb_idx, face_area, normal)
      end do
      call tidy_velocity_fields(u, v)
    end do
  end do

  call fin()

  contains

  !> @brief For a given cell and neighbour computes the global cell and neighbour indices, corresponding face
  !> area, and normal
  !
  !> @param[in] local_idx           - The cell's local index
  !> @param[in] ngb                 - The neighbour we're interested in (range 1-4)
  !> @param[out] self_idx           - The cell's global index
  !> @param[out] ngb_idx            - The neighbour's global index
  !> @param[out] face_area  - The surface area of the face between the cell and its neighbour
  !> @param[out] normal             - The face normal between the cell and its neighbour
  subroutine set_cell_indices(local_idx, ngb, self_idx, ngb_idx, face_area, normal)
    integer(accs_int), intent(in) :: local_idx
    integer(accs_int), intent(in) :: ngb
    integer(accs_int), intent(out) :: self_idx
    integer(accs_int), intent(out) :: ngb_idx
    real(accs_real), intent(out) :: face_area
    real(accs_real), intent(out), dimension(ndim) :: normal

    type(cell_locator) :: self_loc
    type(neighbour_locator) :: ngb_loc
    type(face_locator) :: face_loc
    
    call set_cell_location(self_loc, square_mesh, local_idx)
    call get_global_index(self_loc, self_idx)
  
    call set_neighbour_location(ngb_loc, self_loc, ngb)
    call get_global_index(ngb_loc, ngb_idx)

    call set_face_location(face_loc, square_mesh, local_idx, ngb)
    call get_face_area(face_loc, face_area)

    call get_face_normal(face_loc, normal)
  end subroutine set_cell_indices

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

  !> @brief Checks whether advection coefficient is correct for given velocity fields, cell and neighbour
  !
  !> @param[in] u, v              - The velocity fields
  !> @param[in] self_idx          - The given cell's global index
  !> @param[in] row               - The given cell's row in spatial mesh
  !> @param[in] col               - The given cell's column in spatial mesh
  !> @param[in] ngb_idx           - The neighbour's global index
  !> @param[in] face_area - The surface area of the face between the cell and neighbour
  !> @param[in] normal            - The normal to the face between the cell and neighbour
  subroutine run_advection_coeff_test(u, v, self_idx, row, col, ngb_idx, face_area, normal)
    class(field), intent(in) :: u, v
    integer(accs_int), intent(in) :: self_idx
    integer(accs_int), intent(in) :: row, col
    integer(accs_int), intent(in) :: ngb_idx
    real(accs_real), intent(in) :: face_area
    real(accs_real), dimension(ndim), intent(in) :: normal

    real(accs_real) :: coeff
    real(accs_real) :: expected_coeff

    select type(u)
      type is(central_field)
        select type(v)
          type is(central_field)
            call calc_advection_coeff(self_idx, ngb_idx, face_area, cps, u, v, 0_accs_int, coeff)
          class default
            write(message, *) "FAIL: incorrect velocity field discretisation"
            call stop_test(message)
        end select
      type is(upwind_field)
        select type(v)
          type is(upwind_field)
            call calc_advection_coeff(self_idx, ngb_idx, face_area, cps, u, v, 0_accs_int, coeff)
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
        expected_coeff = min(-face_area*(u%val(col,row)*normal(1) + v%val(col,row)*normal(2)), 0.0_accs_real)
      type is(central_field)
        expected_coeff = -face_area*0.5_accs_real*(u%val(col,row)*normal(1) + v%val(col,row)*normal(2))
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

end program test_advection_coeff