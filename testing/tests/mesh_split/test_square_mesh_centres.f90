!> @brief Test the cell/face centres of a square mesh.
!
!> @description The cell/face centres of a mesh should all fall within the meshed domain, for a
!!              square mesh \f$x\in[0,1]^d\f$.
program test_square_mesh_centres

  use testing_lib

  use constants, only : ndim
  use mesh_utils, only : build_square_mesh, centre

  implicit none
  
  type (mesh) :: square_mesh

  real(accs_real) :: l
  integer(accs_int) :: n

  integer(accs_int) :: i
  integer(accs_int) :: j
  
  type(cell_locator) :: cell_location
  real(accs_real), dimension(ndim) :: cc
  type(face_locator) :: face_location
  real(accs_real), dimension(ndim) :: fc
  
  call init()
  
  do n = 1, 100
    l = parallel_random(par_env)
    square_mesh = build_square_mesh(n, l, par_env)

    do i = 1, square_mesh%nlocal
      call set_cell_location(cell_location, square_mesh, i)
      call centre(cell_location, cc)
      associate(x => cc(1), y => cc(2))
        if ((x > l) .or. (x < 0_accs_real) &
             .or. (y > l) .or. (y < 0_accs_real)) then
          write (message,*) "FAIL: expected cell centre 0 <= x,y <= ", l, " got ", x, " ", y
          call stop_test(message)
        end if
      end associate

      associate(nnb => square_mesh%nnb(i))
        do j = 1, nnb
          call set_face_location(face_location, square_mesh, i, j)
          call centre(face_location, fc)
          associate(x => fc(1), y => fc(2))
            if ((x > (l + eps)) .or. (x < (0.0_accs_real - eps)) &
                 .or. (y > (l + eps)) .or. (y < (0.0_accs_real - eps))) then
              write(message,*) "FAIL: expected face centre 0 <= x,y <= ", l, " got ", x, " ", y
              call stop_test(message)
            end if
          end associate
        end do
      end associate
    end do
    
  end do

  call fin()

end program test_square_mesh_centres
