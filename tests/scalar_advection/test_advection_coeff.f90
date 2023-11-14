!v Test that advection coefficients are calculated correctly
!
!  Computes the advection coefficients for two flow directions
!  (in +x, +y directions) for central and upwind differencing and compares to
!  known values
program test_advection_coeff

  use testing_lib
  use constants, only: ndim, insert_mode
  use types, only: field, upwind_field, central_field, cell_locator, face_locator, neighbour_locator
  use mesh_utils, only: build_square_mesh
  use vec, only: create_vector, get_vector_data, restore_vector_data
  use fv, only: calc_advection_coeff, calc_cell_coords
  use meshing, only: create_cell_locator, create_face_locator, create_neighbour_locator, &
                     get_global_index, get_local_index, get_face_area, get_face_normal, &
                     get_local_num_cells
  use meshing, only: set_mesh_object, nullify_mesh_object
  use utils, only: update, initialise, &
                   set_size, set_row, set_entry, set_values
  use petsctypes, only: vector_petsc

  implicit none
  type(vector_spec) :: vec_properties
  class(field), allocatable :: scalar
  class(field), allocatable :: u, v
  integer(ccs_int), parameter :: cps = 50
  integer(ccs_int) :: index_p, index_nb, index_test
  integer(ccs_int) :: nb
  integer(ccs_int) :: direction, discretisation
  integer, parameter :: x_dir = 1, y_dir = 2
  integer, parameter :: central = -1, upwind = -2
  real(ccs_real) :: face_area
  real(ccs_real), dimension(ndim) :: normal
  real(ccs_real), dimension(:), pointer :: u_data, v_data

  integer(ccs_int) :: local_num_cells
  type(face_locator) :: loc_f

  call init()

  mesh = build_square_mesh(par_env, shared_env, cps, 1.0_ccs_real)
  call set_mesh_object(mesh)

  call get_local_num_cells(local_num_cells)

  index_test = int(0.5 * local_num_cells + 2, ccs_int)
  do direction = x_dir, y_dir
    do discretisation = upwind, central
      if (discretisation == central) then
        allocate (central_field :: scalar)
        allocate (central_field :: u)
        allocate (central_field :: v)
      else if (discretisation == upwind) then
        allocate (upwind_field :: scalar)
        allocate (upwind_field :: u)
        allocate (upwind_field :: v)
      else
        write (message, *) 'Invalid discretisation type selected'
        call stop_test(message)
      end if

      call initialise(vec_properties)
      call set_size(par_env, mesh, vec_properties)
      call create_vector(vec_properties, scalar%values)
      call create_vector(vec_properties, u%values)
      call create_vector(vec_properties, v%values)

      call set_velocity_fields(direction, u, v)

      associate (u_vec => u%values, v_vec => v%values)
        call get_vector_data(u_vec, u_data)
        call get_vector_data(v_vec, v_data)

        do nb = 1, 4
          call get_cell_parameters(index_test, nb, index_p, index_nb, face_area, normal, loc_f)
          call run_advection_coeff_test(scalar, u_data, v_data, index_p, index_nb, face_area, normal, loc_f)
        end do

        call restore_vector_data(u_vec, u_data)
        call restore_vector_data(v_vec, v_data)
      end associate

      call tidy_velocity_fields(scalar, u, v)
    end do
  end do

  call nullify_mesh_object()

  call fin()

contains

  !v For a given cell and neighbour computes the local cell and neighbour indices, corresponding face
  !  area, and normal
  subroutine get_cell_parameters(index, nb, index_p, index_nb, face_area, normal, loc_f)
    integer(ccs_int), intent(in) :: index                  !< The cell's local index
    integer(ccs_int), intent(in) :: nb                     !< The neighbour we're interested in (range 1-4)
    integer(ccs_int), intent(out) :: index_p               !< The cell's local index
    integer(ccs_int), intent(out) :: index_nb              !< The neighbour's local index
    real(ccs_real), intent(out) :: face_area               !< The surface area of the face between the cell and its neighbour
    real(ccs_real), intent(out), dimension(ndim) :: normal !< The face normal between the cell and its neighbour
    type(face_locator), intent(out) :: loc_f               !< Face locator

    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb

    call create_cell_locator(index, loc_p)
    call get_local_index(loc_p, index_p)

    call create_neighbour_locator(loc_p, nb, loc_nb)
    call get_local_index(loc_nb, index_nb)

    call create_face_locator(index, nb, loc_f)
    call get_face_area(loc_f, face_area)

    call get_face_normal(loc_f, normal)
  end subroutine get_cell_parameters

  !v Sets the velocity field in the desired direction and discretisation
  subroutine set_velocity_fields(direction, u, v)

    use vec, only: create_vector_values
    use utils, only: set_mode

    integer(ccs_int), intent(in) :: direction        !< Integer indicating the direction of the velocity field
    class(field), intent(inout), allocatable :: u, v !< The velocity fields in x and y directions
    type(cell_locator) :: loc_p
    type(vector_values) :: u_vals, v_vals
    integer(ccs_int) :: index_p, global_index_p
    real(ccs_real) :: u_val, v_val
    integer(ccs_int) :: n_local

    call get_local_num_cells(n_local)

    call create_vector_values(n_local, u_vals)
    call create_vector_values(n_local, v_vals)
    call set_mode(insert_mode, u_vals)
    call set_mode(insert_mode, v_vals)

    ! Set IC velocity fields
    do index_p = 1, n_local
      call create_cell_locator(index_p, loc_p)
      call get_global_index(loc_p, global_index_p)

      if (direction == x_dir) then
        u_val = 1.0_ccs_real
        v_val = 0.0_ccs_real
      else if (direction == y_dir) then
        u_val = 0.0_ccs_real
        v_val = 1.0_ccs_real
      end if

      call set_row(global_index_p, u_vals)
      call set_entry(u_val, u_vals)

      call set_row(global_index_p, v_vals)
      call set_entry(v_val, v_vals)
    end do

    call set_values(u_vals, u%values)
    call set_values(v_vals, v%values)

    call update(u%values)
    call update(v%values)

    deallocate (u_vals%global_indices)
    deallocate (v_vals%global_indices)
    deallocate (u_vals%values)
    deallocate (v_vals%values)

  end subroutine set_velocity_fields

  !v Deallocates the velocity fields
  subroutine tidy_velocity_fields(scalar, u, v)
    class(field), allocatable :: scalar !< The scalar field structure
    class(field), allocatable :: u, v   !< The velocity fields to deallocate

    deallocate (scalar)
    deallocate (u)
    deallocate (v)
  end subroutine tidy_velocity_fields

  !v Checks whether advection coefficient is correct for given velocity fields, cell and neighbour
  subroutine run_advection_coeff_test(phi, u, v, index_p, index_nb, face_area, face_normal, loc_f)
    class(field), intent(in) :: phi                            !< The scalar field structure
    real(ccs_real), dimension(:), intent(in) :: u, v           !< Arrays containing the velocity fields
    integer(ccs_int), intent(in) :: index_p                    !< The given cell's local index
    integer(ccs_int), intent(in) :: index_nb                   !< The neighbour's local index
    real(ccs_real), intent(in) :: face_area                    !< The surface area of the face between the cell and neighbour
    real(ccs_real), dimension(ndim), intent(in) :: face_normal !< The normal to the face between the cell and neighbour
    type(face_locator), intent(in) :: loc_f                    !< face locator
    
    real(ccs_real) :: coeffaP, coeffaF,coeff
    real(ccs_real) :: mf
    real(ccs_real) :: expected_coeff

    ! Compute mass flux
    mf = 0.5_ccs_real * (u(index_p) + u(index_nb)) * face_normal(1) &
         + 0.5_ccs_real * (v(index_p) + v(index_nb)) * face_normal(2)

    select type (phi)
    type is (central_field)
      call calc_advection_coeff(phi, loc_f, mf, 0_ccs_int, coeffaP, coeffaF)
    type is (upwind_field)
      call calc_advection_coeff(phi, loc_f, mf, 0_ccs_int, coeffaP, coeffaF)
    class default
      write (message, *) "FAIL: incorrect velocity field discretisation"
      call stop_test(message)
    end select
    !coeff = coeff * mf * face_area
    coeff = coeffaF* mf * face_area

    select type (phi)
    type is (upwind_field)
      expected_coeff = min(mf * face_area, 0.0_ccs_real)
    type is (central_field)
      expected_coeff = 0.5_ccs_real * (mf * face_area)
    class default
      write (message, *) "FAIL: incorrect velocity field discretisation"
      call stop_test(message)
    end select

    select type (phi)
    type is (upwind_field)
      if (abs(coeff - expected_coeff) .ge. eps) then
        write (message, *) "FAIL: incorrect upwind advection coefficient computed. Expected ", expected_coeff, " computed ", &
          coeff, " normal ", normal, " u ", u(index_p), " v ", v(index_p)
        call stop_test(message)
      end if
    type is (central_field)
      if (abs(coeff - expected_coeff) .ge. eps) then
        write (message, *) "FAIL: incorrect central advection coefficient computed. Expected ", expected_coeff, " computed ", &
          coeff, " normal ", normal, " u ", u(index_p), " v ", v(index_p)
        call stop_test(message)
      end if
    class default
      write (message, *) "FAIL: incorrect velocity field discretisation"
      call stop_test(message)
    end select

  end subroutine run_advection_coeff_test

end program test_advection_coeff
