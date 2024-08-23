!> @brief Test finding the faces on a specific boundary using the boundary name.
program find_named_faces

  use testing_lib
  use meshing, only: find_entities, set_mesh_object, get_face_normal
  use mesh_utils, only: build_square_mesh
  
  implicit none

  integer, parameter :: cps = 4
  real(ccs_real), parameter :: l = 1.0
  character(len=128), dimension(4) :: bnd_names
  
  type(face_locator), dimension(:), allocatable :: faces

  integer :: i, f

  real(ccs_real), dimension(3) :: n
  
  call init()

  ! Check we are running on one process, this just makes the counting easier rather than anything
  ! else.
  if (par_env%num_procs > 1) then
    call stop_test("This test should be run in serial")
  end if

  ! Construct a boundary name list as though the user had provided one.
  ! XXX: It is the user's responsibility to ensure that this ordering matches the boundary IDs in
  !      the mesh.
  bnd_names(1) = "alpha" ! Left
  bnd_names(2) = "beta"  ! Right
  bnd_names(3) = "gamma" ! Bottom
  bnd_names(4) = "delta" ! Top
  
  mesh = build_square_mesh(par_env, shared_env, cps, l, bnd_names)
  call set_mesh_object(mesh)
  
  do i = 1, size(bnd_names)
    call find_entities(mesh, bnd_names(i), faces)
    if (size(faces) /= cps) then
      call stop_test("Incorrect boundary face count")
    end if

    do f = 1, size(faces)
      call get_face_normal(faces(f), n)

      print *, i, f, n
      if (i == 1) then
        ! Left
        if (n(1) /= -1.0_ccs_real) then
          call stop_test("Expected a left boundary face")
        end if
      else if (i == 2) then
        ! Right
        if (n(1) /= 1.0_ccs_real) then
          call stop_test("Expected a right boundary face")
        end if
      else if (i == 3) then
        ! Bottom
        if (n(2) /= -1.0_ccs_real) then
          call stop_test("Expected a bottom boundary face")
        end if
      else
        ! Top
        if (n(2) /= 1.0_ccs_real) then
          call stop_test("Expected a top boundary face")
        end if
      end if
    end do
  end do
  
  call fin()
  
end program find_named_faces
