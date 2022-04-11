!> @brief Test the square mesh generator creates a closed mesh.
!
!> @description A valid mesh should be "closed" that is the surface integral should be zero. The
!!              same is true of mesh cells.
program test_mesh_square_mesh_closed

  use testing_lib

  use constants

  use meshing, only : set_face_location, get_face_normal, get_face_area
  use mesh_utils, only : build_square_mesh

  implicit none
  
  type(mesh), target :: square_mesh
  type(face_locator) :: face_location

  integer(ccs_int) :: n
  real(ccs_real) :: l

  real(ccs_real), dimension(ndim) :: S
  real(ccs_real), dimension(ndim) :: norm
  real(ccs_real) :: A

  integer(ccs_int) :: i, j

  call init()
  
  do n = 1, 100 ! XXX: Should use some named constant, not just "100"
    l = parallel_random(par_env)
    square_mesh = build_square_mesh(par_env, n, l)

    ! Loop over cells
    do i = 1, square_mesh%nlocal
      S(:) = 0.0_ccs_real

      ! Loop over neighbours/faces
      do j = 1, square_mesh%nnb(i)

        call set_face_location(square_mesh, i, j, face_location)
        call get_face_area(face_location, A)
        call get_face_normal(face_location, norm)
        S(:) = S(:) + norm(:) * A
      end do

      ! Loop over axes
      do j = 1, ndim
        print *, S(j), j
        if (abs(S(j) - 0.0_ccs_real) > eps) then
          write(message, *) "FAIL: expected", 0.0_ccs_real, " got ", S(j)
          call stop_test(message)
        end if
      end do
    end do
  end do

  call fin()
  
end program test_mesh_square_mesh_closed
