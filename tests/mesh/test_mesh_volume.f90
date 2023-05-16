!> @brief Test the square mesh generator creates a correctly-sized mesh.
!
!> @description A "cubic" domain of side L should result in a mesh of volume L^3, this can be
!> verified by summing the volumes of all cells.
program test_mesh_volume

  use testing_lib
  use meshing, only: create_cell_locator, get_volume, get_local_num_cells
  use mesh_utils, only: build_mesh

  implicit none

  type(ccs_mesh) :: mesh

  integer(ccs_int) :: nx, ny, nz
  real(ccs_real) :: l
  real(ccs_real) :: vol
  real(ccs_real) :: vol_global
  real(ccs_real) :: expected_vol

  integer(ccs_int) :: local_num_cells
  integer(ccs_int) :: i
  type(cell_locator) :: loc_p
  real(ccs_real) :: V

  integer(ccs_int) :: nneg_vol
  integer(ccs_int) :: nneg_vol_global

  real(ccs_real) :: CV

  call init()

  nx = 2
  ny = 2
  nz = 2

  l = parallel_random(par_env)
  mesh = build_mesh(par_env, nx, ny, nz, l)
  expected_vol = l**3 ! XXX: Currently the mesh is a hard-coded 3D cube...

  CV = (l / nx) * (l / ny) * (l / nz)

  vol = 0.0_ccs_real
  nneg_vol = 0
  call get_local_num_cells(mesh, local_num_cells)
  do i = 1, local_num_cells
    call create_cell_locator(mesh, i, loc_p)
    call get_volume(loc_p, V)
    if (V <= 0) then
      nneg_vol = nneg_vol + 1
    end if

    if (abs(V - CV) > 1.0e-8) then
      write (message, *) "FAIL: expected cell volume ", CV, " got ", V
      call stop_test(message)
    end if

    vol = vol + V
  end do

  select type (par_env)
  type is (parallel_environment_mpi)
    call MPI_Allreduce(vol, vol_global, 1, real_type, MPI_SUM, par_env%comm, ierr)
    call MPI_Allreduce(nneg_vol, nneg_vol_global, 1, MPI_INT, MPI_SUM, par_env%comm, ierr)
  class default
    write (message, *) "ERROR: Unknown parallel environment!"
    call stop_test(message)
  end select

  ! XXX: This would be a good candidate for a testing library
  if (abs(expected_vol - vol_global) > 1.0e-8) then
    print *, mesh%geo%h, l / nx !TODO: not sure if this should be put inside message
    write (message, *) "FAIL: expected ", expected_vol, " got ", vol_global
    call stop_test(message)
  end if

  if (nneg_vol_global /= 0) then
    write (message, *) "FAIL: encountered negative cell volume!"
    call stop_test(message)
  end if

  call fin()

end program test_mesh_volume
