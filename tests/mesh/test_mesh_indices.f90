!> @brief Test the indexing of cells
program test_mesh_indices

  use testing_lib

  use meshing, only: create_cell_locator, get_global_index, get_local_num_cells, get_global_num_cells
  use mesh_utils, only: build_mesh

  implicit none

  type(ccs_mesh) :: mesh

  real(ccs_real) :: l
  integer(ccs_int) :: n, nx, ny, nz

  integer(ccs_int) :: nlocal
  integer(ccs_int) :: nglobal
  integer(ccs_int) :: i

  type(cell_locator) :: loc_p
  integer(ccs_int) :: global_index

  integer(ccs_int), dimension(5) :: m = (/2, 4, 8, 16, 20/)
  integer(ccs_int) :: mctr

  call init()

  ! XXX: use smaller size than 2D test - 20^3 ~= 100^2
  do mctr = 1, size(m)
    n = m(mctr)

    nx = n
    ny = n
    nz = n

    l = parallel_random(par_env)
    mesh = build_mesh(par_env, nx, ny, nz, l)

    call get_local_num_cells(mesh, nlocal)
    call get_global_num_cells(mesh, nglobal)
    do i = 1, nlocal
      call create_cell_locator(mesh, i, loc_p)
      call get_global_index(loc_p, global_index)
      if ((global_index < 1) .or. (global_index > nglobal)) then
        if (global_index /= -1) then
          write (message, *) "FAIL: expected global index 1 <= idx <= ", nglobal, " got ", global_index
          call stop_test(message)
        end if
        exit
      end if
    end do
  end do

  call fin()

end program test_mesh_indices
