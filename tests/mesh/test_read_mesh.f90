!> @brief Read mesh file and perform some checks on the resulting mesh object:
! - Check that the sum of cell volume is the 'theoretical' mesh volume (pi**3 for the TGV mesh)
! - Check that for each cell its neighbour connectivity is reciprocal
! - Check that for each cell normal summation add up to 0
program test_read_mesh

  use testing_lib
  use core, only: ccs_options, read_input_mesh
  use ccs_base, only: bnd_names_default
  use ccs_base, only: mesh
  use mesh_utils, only: read_mesh, test_mesh_internal_neighbours
  use meshing, only: set_mesh_object, nullify_mesh_object
  use meshing, only: create_cell_locator, create_neighbour_locator, count_neighbours, &
                     get_boundary_status, get_local_num_cells, get_local_index, create_face_locator

  implicit none

  type(ccs_options):: run_options

  
  call init()

  run_options%paths%case_name = "MESH_FILES/cube_tet_128k"
  run_options%mesh%bnd_names = bnd_names_default
  
  call read_mesh(par_env, shared_env, run_options, mesh)
  call set_mesh_object(mesh)

  ! Check mesh volume
  call check_mesh_volume()

  ! Check mesh neighbour connectivity
  call check_neighbours()

  ! Check normal integration is 0
  call check_normals()
 
  call nullify_mesh_object()
  call fin()

  contains

  !v Checks that for each cell, the sum of its face normals is 0
  subroutine check_normals()

    use meshing, only: get_face_normal

    integer(ccs_int) :: nnb
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: i, j
    logical :: is_boundary
    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb
    type(face_locator) :: loc_f 
    real(ccs_real), dimension(3) :: face_normal
    real(ccs_real), dimension(3) :: face_normal_sum

    call get_local_num_cells(local_num_cells)
    do i = 1, local_num_cells

      call create_cell_locator(i, loc_p)
      call count_neighbours(loc_p, nnb)
      face_normal_sum(:) = 0.0_ccs_real

      ! Loop over neighbours
      do j = 1, nnb
        call create_face_locator(i, j, loc_f)
        call create_neighbour_locator(loc_p, j, loc_nb)
        call get_face_normal(loc_f, face_normal)
        call get_boundary_status(loc_nb, is_boundary)

        if (is_boundary) then
          ! Boundary faces are inwards facing
          face_normal_sum = face_normal_sum - face_normal
        else
          face_normal_sum = face_normal_sum + face_normal
        end if
      end do

      call assert_eq(norm2(face_normal_sum), 0.0_ccs_real, "Sum of face normals not zero")

    end do

  end subroutine


  !v Check that each neighbour are linked together. Makes sure topology is consistent
  subroutine check_neighbours()

    integer(ccs_int) :: nnb
    integer(ccs_int) :: index_nb
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: i, j
    logical :: is_boundary
    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb

    call get_local_num_cells(local_num_cells)
    do i = 1, local_num_cells

      call create_cell_locator(i, loc_p)
      call count_neighbours(loc_p, nnb)

      ! Loop over neighbours
      do j = 1, nnb
        call create_neighbour_locator(loc_p, j, loc_nb)

        call get_local_index(loc_nb, index_nb)
        call assert_neq(index_nb, 0, "All neighbours should be filled!")
        
        call get_boundary_status(loc_nb, is_boundary)
        if (.not. is_boundary) then
          call test_mesh_internal_neighbours(loc_nb)
        end if
      end do

    end do

  end subroutine 


  !v Check that the sum of cell volume adds up to the mesh overall volume
  subroutine check_mesh_volume()
  use meshing, only: create_cell_locator, get_volume, get_local_num_cells

    real(ccs_real):: vol, V, vol_global
    integer(ccs_int):: i, local_num_cells, nneg_vol, nneg_vol_global
    type(cell_locator):: loc_p

    integer:: ierr
    
    vol = 0.0_ccs_real
    nneg_vol = 0
    call get_local_num_cells(local_num_cells)
    do i = 1, local_num_cells
        call create_cell_locator(i, loc_p)
        call get_volume(loc_p, V)
        if (V <= 0) then
            nneg_vol = nneg_vol+1
        end if

        vol = vol+V
    end do

    select type (par_env)
    type is (parallel_environment_mpi)
        call MPI_Allreduce(vol, vol_global, 1, real_type, MPI_SUM, par_env%comm, ierr)
        call MPI_Allreduce(nneg_vol, nneg_vol_global, 1, MPI_INT, MPI_SUM, par_env%comm, ierr)
    class default
        write (message, *) "ERROR: Unknown parallel environment!"
        call stop_test(message)
    end select

    ! 31.006.... ~= pi**3
    call assert_eq(vol_global, 31.006280191406404d0, "Mesh volume not matching")
    call assert_eq(nneg_vol_global, 0, "Cells with negative volume detected")

  end subroutine

end program test_read_mesh
