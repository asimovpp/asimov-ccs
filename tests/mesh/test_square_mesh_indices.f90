!> @brief Test the indexing of cells
program test_square_mesh_indices

  use testing_lib

  use meshing, only : set_cell_location, get_global_index
  use mesh_utils, only : build_square_mesh

  implicit none

  type(ccs_mesh) :: mesh

  real(ccs_real) :: l
  integer(ccs_int) :: n

  integer(ccs_int) :: i

  type(cell_locator) :: cell_location
  integer(ccs_int) :: idxg

  call init()
  
  do n = 1, 100 ! TODO: Investigate how we can replicate nmax=100 across multiple test programs
    l = parallel_random(par_env)
    mesh = build_square_mesh(par_env, n, l)

    associate(nlocal => mesh%nlocal, &
         nglobal => mesh%nglobal)
      do i = 1, nlocal
        call set_cell_location(mesh, i, cell_location)
        call get_global_index(cell_location, idxg)
        if ((idxg < 1) .or. (idxg > nglobal)) then
          if (idxg /= -1) then
            write(message, *) "FAIL: expected global index 1 <= idx <= ", nglobal, " got ", idxg
            call stop_test(message)
          end if
          exit
        end if
      end do
    end associate
  end do

  call fin()
  
end program test_square_mesh_indices
