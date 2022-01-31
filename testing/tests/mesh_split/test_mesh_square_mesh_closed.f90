!> @brief Test the square mesh generator creates a closed mesh.
!
!> @description A valid mesh should be "closed" that is the surface integral should be zero. The
!!              same is true of mesh cells.
program test_mesh_square_mesh_closed

  use testing_lib

  use constants

  use meshing, only : set_face_location
  use mesh_utils, only : build_square_mesh, face_normal, face_area

  implicit none
  
  type(mesh), target :: square_mesh
  type(face_locator) :: face_location

  integer(accs_int) :: n
  real(accs_real) :: l

  real(accs_real), dimension(ndim) :: S
  real(accs_real), dimension(ndim) :: norm
  real(accs_real) :: A

  integer(accs_int) :: i, j

  call init()
  
  do n = 1, 100 ! XXX: Should use some named constant, not just "100"
    l = parallel_random(par_env)
    square_mesh = build_square_mesh(n, l, par_env)

    ! Loop over cells
    do i = 1, square_mesh%nlocal
      S(:) = 0.0_accs_real

      ! Loop over neighbours/faces
      do j = 1, square_mesh%nnb(i)

        call set_face_location(face_location, square_mesh, i, j)
        call face_area(face_location, A)
        call face_normal(face_location, norm)
        S(:) = S(:) + norm(:) * A
      end do

      ! Loop over axes
      do j = 1, ndim
        print *, S(j), j
        if (abs(S(j) - 0.0_accs_real) > eps) then
          write(message, *) "FAIL: expected", 0.0_accs_real, " got ", S(j)
          call stop_test(message)
        end if
      end do
    end do
  end do

  call fin()
  
end program test_mesh_square_mesh_closed
