!> @brief Test the indexing of cells
program test_mesh_indices

  use testing_lib

  use ccs_base, only: bnd_names_default
  use core
  
  use meshing, only: create_cell_locator, get_global_index, get_local_num_cells, get_global_num_cells
  use meshing, only: set_mesh_object, nullify_mesh_object
  use mesh_utils, only: build_mesh

  implicit none

  real(ccs_real) :: l
  integer(ccs_int) :: n

  integer(ccs_int) :: nlocal
  integer(ccs_int) :: nglobal
  integer(ccs_int) :: i

  type(cell_locator) :: loc_p
  integer(ccs_int) :: global_index

  integer(ccs_int), dimension(5) :: m = (/4, 8, 12 ,16, 20/)
  integer(ccs_int) :: mctr

  type(ccs_options) :: run_options

  call init()

  ! XXX: use smaller size than 2D test - 20^3 ~= 100^2
  do mctr = 2, size(m)
    n = m(mctr)

    l = parallel_random(par_env)
    run_options%mesh%bnd_names = bnd_names_default
    run_options%mesh%cps = n
    run_options%mesh%domain_size = l
    mesh = build_mesh(par_env, shared_env, run_options)
    call set_mesh_object(mesh)

    call get_local_num_cells(nlocal)
    call get_global_num_cells(nglobal)
    do i = 1, nlocal
      call create_cell_locator(i, loc_p)
      call get_global_index(loc_p, global_index)
      if ((global_index < 1) .or. (global_index > nglobal)) then
        if (global_index /= -1) then
          write (message, *) "FAIL: expected global index 1 <= idx <= ", nglobal, " got ", global_index
          call stop_test(message)
        end if
        exit
      end if
    end do

    call nullify_mesh_object()
  end do

  call fin()

end program test_mesh_indices
