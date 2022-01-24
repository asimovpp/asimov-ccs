!> @brief Test the indexing of cells
program test_square_mesh_indices

  use testing_lib
  
  use mesh_utils, only : build_square_mesh, global_index

  implicit none

  type (mesh) :: square_mesh

  real(accs_real) :: l
  integer(accs_int) :: n

  integer(accs_int) :: i

  type(cell_locator) :: cell_location
  integer(accs_int) :: idxg

  do n = 1, 100 ! TODO: Investigate how we can replicate nmax=100 across multiple test programs
    l = parallel_random(par_env)
    square_mesh = build_square_mesh(n, l, par_env)

    associate(nlocal => square_mesh%nlocal, &
         nglobal => square_mesh%n)
      do i = 1, nlocal
        call set_cell_location(cell_location, square_mesh, i)
        call global_index(cell_location, idxg)
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

end program test_square_mesh_indices
