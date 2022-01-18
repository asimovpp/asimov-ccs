!> @brief Test the cell centres of a square mesh.
!
!> @description The cell centres of a mesh should all fall within the meshed domain, for a square
!!              mesh \f$x\in[0,1]^d\f$.
program test_square_mesh_cell_centres

  use testing_lib
  use mesh_utils, only : build_square_mesh

  type (mesh) :: square_mesh

  real(accs_real) :: l
  integer(accs_int) :: n

  integer(accs_int) :: i
  
  call init()
  
  do n = 1, 100
    l = parallel_random(par_env)
    square_mesh = build_square_mesh(n, l, par_env)

    do i = 1, square_mesh%nlocal
      associate(x => square_mesh%xc(1, i), &
           y => square_mesh%xc(2, i))
        if ((x > l) .or. (x < 0_accs_real) &
             .or. (y > l) .or. (y < 0_accs_real)) then
          write (message,*) "FAIL: expected 0 <= x,y <= ", l, " got ", x, " ", y
          call stop_test(message)
        end if
      end associate
    end do
    
  end do

  call fin()

end program test_square_mesh_cell_centres
