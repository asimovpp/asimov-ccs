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
                     get_global_index, get_local_index, get_face_area, get_face_normal
  use utils, only : update, initialise, &
                set_size, pack_entries, set_values
  use petsctypes, only: vector_petsc

  type(ccs_mesh) :: square_mesh
  type(vector_spec) :: vec_sizes
  class(field), allocatable :: scalar
  class(field), allocatable :: u, v
  integer(ccs_int), parameter :: cps = 50
  integer(ccs_int) :: self_idx, index_nb, local_idx
  integer(ccs_int) :: ngb
  integer(ccs_int) :: direction, discretisation
  integer, parameter :: x_dir = 1, y_dir = 2
  integer, parameter :: central = -1, upwind = -2
  real(ccs_real) :: face_area
  real(ccs_real), dimension(ndim) :: normal
  real(ccs_real), dimension(:), pointer :: u_data, v_data
  
  call init()

  square_mesh = build_square_mesh(par_env, cps, 1.0_ccs_real)

  local_idx = int(0.5*square_mesh%nlocal + 2, ccs_int)
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
      call set_size(par_env, square_mesh, vec_sizes)
      call create_vector(vec_sizes, scalar%values)
      call create_vector(vec_sizes, u%values)
      call create_vector(vec_sizes, v%values)

      call set_velocity_fields(square_mesh, direction, u, v)

      associate (u_vec => u%values, v_vec => v%values)
        call get_vector_data(u_vec, u_data)
        call get_vector_data(v_vec, v_data)

        do ngb = 1, 4
          call get_cell_parameters(local_idx, ngb, self_idx, index_nb, face_area, normal)
          call run_advection_coeff_test(scalar, u_data, v_data, self_idx, index_nb, face_area, normal)
        end do
      
        call restore_vector_data(u_vec, u_data)
        call restore_vector_data(v_vec, v_data)
      end associate

      call tidy_velocity_fields(scalar, u, v)
    end do
  end do

  call fin()

  contains

  !> @brief For a given cell and neighbour computes the local cell and neighbour indices, corresponding face
  !> area, and normal
  !
  !> @param[in] local_idx           - The cell's local index
  !> @param[in] ngb                 - The neighbour we're interested in (range 1-4)
  !> @param[out] self_idx           - The cell's local index
  !> @param[out] index_nb            - The neighbour's local index
  !> @param[out] face_area  - The surface area of the face between the cell and its neighbour
  !> @param[out] normal             - The face normal between the cell and its neighbour
  subroutine get_cell_parameters(local_idx, ngb, self_idx, index_nb, face_area, normal)
    integer(ccs_int), intent(in) :: local_idx
    integer(ccs_int), intent(in) :: ngb
    integer(ccs_int), intent(out) :: self_idx
    integer(ccs_int), intent(out) :: index_nb
    real(ccs_real), intent(out) :: face_area
    real(ccs_real), intent(out), dimension(ndim) :: normal

    type(cell_locator) :: self_loc
    type(neighbour_locator) :: loc_nb
    type(face_locator) :: face_loc
    
    call set_cell_location(square_mesh, local_idx, self_loc)
    call get_local_index(self_loc, self_idx)
  
    call set_neighbour_location(self_loc, ngb, loc_nb)
    call get_local_index(loc_nb, index_nb)

    call set_face_location(square_mesh, local_idx, ngb, face_loc)
    call get_face_area(face_loc, face_area)

    call get_face_normal(face_loc, normal)
  end subroutine get_cell_parameters

  !> @brief Sets the velocity field in the desired direction and discretisation
  !
  !> @param[in] cell_mesh - The mesh structure
  !> @param[in] direction - Integer indicating the direction of the velocity field
  !> @param[out] u, v     - The velocity fields in x and y directions
  subroutine set_velocity_fields(cell_mesh, direction, u, v)
    class(ccs_mesh), intent(in) :: cell_mesh
    integer(ccs_int), intent(in) :: direction
    class(field), intent(inout), allocatable :: u, v
    type(cell_locator) :: self_loc
    type(vector_values) :: u_vals, v_vals
    integer(ccs_int) :: local_idx, self_idx
    real(ccs_real) :: u_val, v_val

    u_vals%setter_mode = insert_mode
    v_vals%setter_mode = insert_mode
    
    associate(n_local => cell_mesh%nlocal)
      allocate(u_vals%indices(n_local))
      allocate(v_vals%indices(n_local))
      allocate(u_vals%values(n_local))
      allocate(v_vals%values(n_local))
      
      ! Set IC velocity fields
      do local_idx = 1, n_local
        call set_cell_location(cell_mesh, local_idx, self_loc)
        call get_global_index(self_loc, self_idx)

        if (direction == x_dir) then
          u_val = 1.0_ccs_real
          v_val = 0.0_ccs_real
        else if (direction == y_dir) then
          u_val = 0.0_ccs_real
          v_val = 1.0_ccs_real
        end if

        call pack_entries(local_idx, self_idx, u_val, u_vals)
        call pack_entries(local_idx, self_idx, v_val, v_vals)
      end do
    end associate
    call set_values(u_vals, u%values)
    call set_values(v_vals, v%values)

    call update(u%values)
    call update(v%values)
    
    deallocate(u_vals%indices)
    deallocate(v_vals%indices)
    deallocate(u_vals%values)
    deallocate(v_vals%values)
    
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
  !> @param[in] self_idx    - The given cell's local index
  !> @param[in] index_nb     - The neighbour's local index
  !> @param[in] face_area   - The surface area of the face between the cell and neighbour
  !> @param[in] face_normal - The normal to the face between the cell and neighbour
  subroutine run_advection_coeff_test(phi, u, v, self_idx, index_nb, face_area, face_normal)
    class(field), intent(in) :: phi
    real(ccs_real), dimension(:), intent(in) :: u, v
    integer(ccs_int), intent(in) :: self_idx
    integer(ccs_int), intent(in) :: index_nb
    real(ccs_real), intent(in) :: face_area
    real(ccs_real), dimension(ndim), intent(in) :: face_normal

    real(ccs_real) :: coeff
    real(ccs_real) :: mf
    real(ccs_real) :: expected_coeff

    !! Compute mass flux
    mf = 0.5_ccs_real * (u(self_idx) + u(index_nb)) * face_normal(1) &
         + 0.5_ccs_real * (v(self_idx) + v(index_nb)) * face_normal(2)
    
    select type(phi)
      type is(central_field)
        call calc_advection_coeff(phi, mf, 0_ccs_int, coeff)
      type is(upwind_field)
        call calc_advection_coeff(phi, mf, 0_ccs_int, coeff)
      class default
        write(message, *) "FAIL: incorrect velocity field discretisation"
        call stop_test(message)
    end select
    coeff = coeff * mf * face_area 

    select type(phi)
      type is(upwind_field)
        expected_coeff = min(mf * face_area, 0.0_ccs_real)
      type is(central_field)
        expected_coeff = 0.5_ccs_real * (mf * face_area)
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
