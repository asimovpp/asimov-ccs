!> @brief Test that advection coefficients are calculated correctly
!
!> @description Computes the advection coefficients for two flow directions
!> (in +x, +y directions) for central and upwind differencing and compares to
!> known values
program test_advection_coeff

  use testing_lib
  use constants, only: ndim, insert_mode
  use types, only: field, upwind_field, central_field, cell_locator, face_locator, neighbour_locator
  use mesh_utils, only : build_square_mesh
  use vec, only : create_vector, get_vector_data, restore_vector_data
  use fv, only: calc_advection_coeff, calc_cell_coords
  use meshing, only: set_cell_location, set_face_location, set_neighbour_location, &
                     get_global_index, get_face_area, get_face_normal
  use utils, only : update, initialise, &
                set_global_size, pack_entries, set_values
  use petsctypes, only: vector_petsc

  type(mesh) :: square_mesh
  type(vector_init_data) :: vec_sizes
  class(field), allocatable :: scalar
  class(field), allocatable :: u, v
  integer(accs_int), parameter :: cps = 50
  integer(accs_int) :: self_idx, ngb_idx, local_idx
  integer(accs_int) :: ngb
  integer(accs_int) :: direction, discretisation
  integer, parameter :: x_dir = 1, y_dir = 2
  integer, parameter :: central = -1, upwind = -2
  real(accs_real) :: face_area
  real(accs_real), dimension(ndim) :: normal
  real(accs_real), dimension(:), pointer :: u_data, v_data
  
  call init()

  square_mesh = build_square_mesh(cps, 1.0_accs_real, par_env)

  local_idx = int(0.5*square_mesh%nlocal + 2, accs_int)
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

      associate (u_vec => u%vec, v_vec => v%vec)
        call get_vector_data(u_vec, u_data)
        call get_vector_data(v_vec, v_data)

        do ngb = 1, 4
          call get_cell_parameters(local_idx, ngb, self_idx, ngb_idx, face_area, normal)
          call run_advection_coeff_test(scalar, u_data, v_data, self_idx, ngb_idx, face_area, normal)
        end do
      
        call restore_vector_data(u_vec, u_data)
        call restore_vector_data(v_vec, v_data)
      end associate

      call tidy_velocity_fields(scalar, u, v)
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
  subroutine get_cell_parameters(local_idx, ngb, self_idx, ngb_idx, face_area, normal)
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
  end subroutine get_cell_parameters

  !> @brief Sets the velocity field in the desired direction and discretisation
  !
  !> @param[in] cell_mesh - The mesh structure
  !> @param[in] direction - Integer indicating the direction of the velocity field
  !> @param[out] u, v     - The velocity fields in x and y directions
  subroutine set_velocity_fields(cell_mesh, direction, u, v)
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

    call update(u%vec)
    call update(v%vec)
    
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

  !> @brief Checks whether advection coefficient is correct for given velocity fields, cell and neighbour
  !
  !> @param[in] scalar      - The scalar field structure
  !> @param[in] u, v        - Arrays containing the velocity fields
  !> @param[in] self_idx    - The given cell's global index
  !> @param[in] ngb_idx     - The neighbour's global index
  !> @param[in] face_area   - The surface area of the face between the cell and neighbour
  !> @param[in] face_normal - The normal to the face between the cell and neighbour
  subroutine run_advection_coeff_test(phi, u, v, self_idx, ngb_idx, face_area, face_normal)
    class(field), intent(in) :: phi
    real(accs_real), dimension(:), intent(in) :: u, v
    integer(accs_int), intent(in) :: self_idx
    integer(accs_int), intent(in) :: ngb_idx
    real(accs_real), intent(in) :: face_area
    real(accs_real), dimension(ndim), intent(in) :: face_normal

    real(accs_real) :: coeff
    real(accs_real) :: mf
    real(accs_real) :: expected_coeff

    !! Compute mass flux
    mf = 0.5_accs_real * (u(self_idx) + u(ngb_idx)) * normal(1) &
         + 0.5_accs_real * (v(self_idx) + v(ngb_idx)) * normal(2)
    
    select type(phi)
      type is(central_field)
        call calc_advection_coeff(phi, ngb_idx, self_idx, face_area, face_normal, u, v, 0_accs_int, coeff)
      type is(upwind_field)
        call calc_advection_coeff(phi, ngb_idx, self_idx, face_area, face_normal, u, v, 0_accs_int, coeff)
      class default
        write(message, *) "FAIL: incorrect velocity field discretisation"
        call stop_test(message)
    end select

    select type(phi)
      type is(upwind_field)
        expected_coeff = min(mf * face_area, 0.0_accs_real)
      type is(central_field)
        expected_coeff = 0.5_accs_real * (mf * face_area)
      class default
        write(message, *) "FAIL: incorrect velocity field discretisation"
        call stop_test(message)
    end select
        
    select type(phi)
      type is(upwind_field)
        if (abs(coeff - expected_coeff) .ge. eps) then
          write(message, *) "FAIL: incorrect upwind advection coefficient computed. Expected ", expected_coeff, " computed ", &
                            coeff, " normal ", normal, " u ", u(self_idx), " v ", v(self_idx)
          call stop_test(message)
        end if
      type is(central_field)
        if (abs(coeff - expected_coeff) .ge. eps) then
          write(message, *) "FAIL: incorrect central advection coefficient computed. Expected ", expected_coeff, " computed ", &
                            coeff, " normal ", normal, " u ", u(self_idx), " v ", v(self_idx)
          call stop_test(message)
        end if
      class default
        write(message, *) "FAIL: incorrect velocity field discretisation"
        call stop_test(message)
    end select

  end subroutine run_advection_coeff_test

end program test_advection_coeff
