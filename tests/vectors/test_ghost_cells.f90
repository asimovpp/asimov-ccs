!> @brief Test setting of ghosts cells
program test_ghost_cells

  use testing_lib
  use constants, only: insert_mode
  use kinds, only: ccs_int
  use types, only: field, upwind_field, central_field, cell_locator, face_locator, neighbour_locator
  use mesh_utils, only: build_square_mesh
  use vec, only: create_vector, update_vector, get_vector_data, restore_vector_data
  use meshing, only: create_neighbour_locator, create_cell_locator, &
                     get_global_index, get_local_index, get_face_area, get_face_normal, &
                     get_local_num_cells, &
                     get_total_num_cells
  use utils, only: update, initialise, &
                   set_size, set_values
  use petsctypes, only: vector_petsc

  implicit none

  type(vector_spec) :: vec_properties
  type(ccs_mesh) :: mesh
  class(ccs_vector), allocatable :: v
  real(ccs_real), dimension(:), pointer :: values

  integer(ccs_int) :: local_num_cells
  integer(ccs_int) :: total_num_cells
  integer(ccs_int) :: i
  integer(ccs_int) :: proc_id
  integer(ccs_int) :: num_procs

  type(cell_locator) :: loc_p
  integer(ccs_int) :: global_index_p

  call init()

  proc_id = par_env%proc_id
  num_procs = par_env%num_procs

  mesh = build_square_mesh(par_env, shared_env, 11, 1.0_ccs_real)
  call get_local_num_cells(mesh, local_num_cells)

  ! Specify vector size based on the mesh
  call initialise(vec_properties)
  call set_size(par_env, mesh, vec_properties)

  ! Create the vector
  call create_vector(vec_properties, v)

  ! Retrieve initial vector values
  call get_vector_data(v, values)

  ! Set vector values to global mesh indices
  do i = 1, local_num_cells
    call create_cell_locator(mesh, i, loc_p)
    call get_global_index(loc_p, global_index_p)
    values(i) = global_index_p
  end do

  ! Restore vector data
  call restore_vector_data(v, values)

  ! Now update the vector (including ghost cells)
  call update(v)

  ! Retrieve the new vector values (including ghost cells)
  call get_vector_data(v, values)

  call get_total_num_cells(mesh, total_num_cells)
  do i = 1, total_num_cells
    call create_cell_locator(mesh, i, loc_p)
    call get_global_index(loc_p, global_index_p)
    if (values(i) /= global_index_p) then
      write (message, *) 'FAIL: wrong vector value. Expected ', global_index_p, ', got ', values(i)
      call stop_test(message)
    end if
  end do

  ! Remove reference to the values array
  call restore_vector_data(v, values)

  call fin()

end program test_ghost_cells
