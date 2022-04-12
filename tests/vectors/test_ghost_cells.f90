!> @brief Test setting of ghosts cells
program test_ghost_cells

  use testing_lib
  use constants, only: insert_mode
  use kinds, only: ccs_int
  use types, only: field, upwind_field, central_field, cell_locator, face_locator, neighbour_locator
  use mesh_utils, only : build_square_mesh
  use vec, only : create_vector, update_vector, get_vec_properties, restore_vec_properties
  use meshing, only: set_neighbour_location, &
                     get_global_index, get_local_index, get_face_area, get_face_normal
  use utils, only : update, initialise, &
                set_size, pack_entries, set_values
  use petsctypes, only: vector_petsc

  implicit none

  type(vector_spec) :: vec_properties
  type(ccs_mesh) :: cell_mesh
  class(ccs_vector), allocatable :: v
  real(ccs_real), dimension(:), pointer :: values

  integer(ccs_int) :: i
  integer(ccs_int) :: proc_id
  integer(ccs_int) :: num_procs

  call init()

  proc_id = par_env%proc_id
  num_procs = par_env%num_procs

  cell_mesh = build_square_mesh(par_env, 11, 1.0_ccs_real)

  ! Specify vector size based on the mesh
  call initialise(vec_properties)
  call set_size(par_env, cell_mesh, vec_properties)

  ! Create the vector
  call create_vector(vec_properties, v)

  ! Retried initial vector values
  call get_vec_properties(v, values)

  ! Set vector values to global mesh indices
  do i = 1, cell_mesh%nlocal
    values(i) = cell_mesh%idx_global(i)
  end do

  ! Now update the vector (including ghost cells)
  call update(v)

  ! Retrieve the new vector values (including ghost cells)
  call get_vec_properties(v, values)

  do i = 1, cell_mesh%ntotal
    if(values(i) /= cell_mesh%idx_global(i)) then
      write(message, *) 'FAIL: wrong vector value. Expected ', cell_mesh%idx_global(i), ', got ', values(i)
      call stop_test(message)
    end if
  end do

  ! Remove reference to the values array
  call restore_vec_properties(v, values)

  call fin()

end program test_ghost_cells
