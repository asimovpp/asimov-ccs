!v Test the indexing of cells
program test_square_mesh_indices

  use testing_lib

  use meshing, only: set_cell_location, get_global_index
  use mesh_utils, only: build_square_mesh

  implicit none

  type(ccs_mesh) :: mesh

  real(ccs_real) :: l
  integer(ccs_int) :: n

  integer(ccs_int) :: i

  type(cell_locator) :: loc_p
  integer(ccs_int) :: global_index

  call init()

  do n = 1, 100 ! TODO: Investigate how we can replicate nmax=100 across multiple test programs
    l = parallel_random(par_env)
    mesh = build_square_mesh(par_env, n, l)

    associate (nlocal => mesh%topo%local_num_cells, &
               nglobal => mesh%topo%global_num_cells)
      do i = 1, nlocal
        call set_cell_location(mesh, i, loc_p)
        call get_global_index(loc_p, global_index)
        if ((global_index < 1) .or. (global_index > nglobal)) then
          if (global_index /= -1) then
            write (message, *) "FAIL: expected global index 1 <= idx <= ", nglobal, " got ", global_index
            call stop_test(message)
          end if
          exit
        end if
      end do
    end associate
  end do

  call fin()

end program test_square_mesh_indices
