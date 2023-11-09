!> @brief Test the cell/face centres of a mesh.
!
!> @description The cell/face centres of a mesh should all fall within the meshed domain, for a
!!              square mesh \f$x\in[0,1]^d\f$.
program test_mesh_centres

  use testing_lib

  use constants, only: ndim
  use meshing, only: create_cell_locator, create_face_locator, create_vert_locator, get_centre, &
                     get_local_num_cells, get_vert_per_cell
  use mesh_utils, only: build_mesh
  use meshing, only: set_mesh_object, nullify_mesh_object

  implicit none


  real(ccs_real) :: l
  integer(ccs_int) :: n, nx, ny, nz

  integer(ccs_int) :: local_num_cells
  integer(ccs_int) :: vert_per_cell
  integer(ccs_int) :: i
  integer(ccs_int) :: j

  type(cell_locator) :: loc_p
  real(ccs_real), dimension(ndim) :: cc
  type(face_locator) :: loc_f
  real(ccs_real), dimension(ndim) :: fc
  type(vert_locator) :: loc_v
  real(ccs_real), dimension(ndim) :: vc

  integer :: dim

  integer(ccs_int), dimension(5) :: m = (/4, 8, 12, 16, 20/)
  integer(ccs_int) :: mctr

  call init()

  ! XXX: use smaller size than 2D test - 20^3 ~= 100^2
  do mctr = 2, size(m)
    n = m(mctr)

    nx = n
    ny = n
    nz = n

    l = parallel_random(par_env)
    mesh = build_mesh(par_env, shared_env, nx, ny, nz, l)
    call set_mesh_object(mesh)

    call get_local_num_cells(local_num_cells)
    do i = 1, local_num_cells
      call create_cell_locator(i, loc_p)
      call get_centre(loc_p, cc)
      associate (x => cc(1), y => cc(2))
        if ((x > l) .or. (x < 0_ccs_real) &
            .or. (y > l) .or. (y < 0_ccs_real)) then
          write (message, *) "FAIL: expected cell centre 0 <= x,y <= ", l, " got ", x, " ", y
          call stop_test(message)
        end if
      end associate

      associate (nnb => mesh%topo%num_nb(i))
        do j = 1, nnb
          call create_face_locator(i, j, loc_f)
          call get_centre(loc_f, fc)
          associate (x => fc(1), y => fc(2))
            if ((x > (l + eps)) .or. (x < (0.0_ccs_real - eps)) &
                .or. (y > (l + eps)) .or. (y < (0.0_ccs_real - eps))) then
              write (message, *) "FAIL: expected face centre 0 <= x,y <= ", l, " got ", x, " ", y
              call stop_test(message)
            end if
          end associate
        end do
      end associate

      call get_vert_per_cell(vert_per_cell)

      do j = 1, vert_per_cell
        call create_vert_locator(i, j, loc_v)
        call get_centre(loc_v, vc)
        do dim = 1, ndim
          if ((vc(dim) > (l + eps)) .or. (vc(dim) < (0.0_ccs_real - eps))) then
            write (message, *) "FAIL: expected vertex centre in range [0, ", l, "] got ", vc
            call stop_test(message)
          end if
        end do
      end do
    end do
  end do

  call nullify_mesh_object()

  call fin()

end program test_mesh_centres
