!> @brief Test setting of ghosts cells
program test_ghost_cells

  use testing_lib
  use constants, only: insert_mode
  use kinds, only: accs_int
  use types, only: field, upwind_field, central_field, cell_locator, face_locator, neighbour_locator
  use mesh_utils, only : build_square_mesh
  use vec, only : create_vector, update_vector, get_vector_data, restore_vector_data
  use meshing, only: set_cell_location, set_face_location, set_neighbour_location, &
                     get_global_index, get_local_index, get_face_area, get_face_normal
  use utils, only : update, initialise, &
                set_global_size, pack_entries, set_values
  use petsctypes, only: vector_petsc

  implicit none

  type(vector_init_data) :: vector_data
  type(mesh) :: cell_mesh
  class(vector), allocatable :: v
  real(accs_real), dimension(:), pointer :: values

  integer(accs_int) :: i
  integer(accs_int) :: proc_id
  integer(accs_int) :: num_procs

  call init()

  proc_id = par_env%proc_id
  num_procs = par_env%num_procs

  cell_mesh = build_square_mesh(11, 1.0_accs_real, par_env)

  ! Specify vector size based on the mesh
  call initialise(vector_data)
  call set_global_size(vector_data, cell_mesh, par_env)

  ! Create the vector
  call create_vector(vector_data, v)

  ! Retried initial vector values
  call get_vector_data(v, values)

  ! Set vector values to global mesh indices
  do i = 1, cell_mesh%nlocal
    values(i) = cell_mesh%idx_global(i)
  end do

  ! Now update the vector (including ghost cells)
  call update(v)

  ! Retrieve the new vector values (including ghost cells)
  call get_vector_data(v, values)

  do i = 1, cell_mesh%ntotal
    if(values(i) /= cell_mesh%idx_global(i)) then
      write(message, *) 'FAIL: wrong vector value. Expected ', cell_mesh%idx_global(i), ', got ', values(i)
      call stop_test(message)
    end if
  end do

  ! Remove reference to the values array
  call restore_vector_data(v, values)

  call fin()

end program test_ghost_cells
