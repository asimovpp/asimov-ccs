!> @brief Read mesh file and check if that has been done properly
program test_read_mesh

  use testing_lib
  use core, only: ccs_options, read_input_mesh
  use ccs_base, only: bnd_names_default
  use ccs_base, only: mesh
  use mesh_utils, only: read_mesh
  use meshing, only: set_mesh_object, nullify_mesh_object

  implicit none

  type(ccs_options):: run_options

  
  call init()

  run_options%paths%case_name = "MESH_FILES/cube_tet_128k"
  run_options%mesh%bnd_names = bnd_names_default
  
  call read_mesh(par_env, shared_env, run_options, mesh)
  call set_mesh_object(mesh)



  ! Check mesh volume
  call check_mesh_volume()


!   call check_mesh_neighbours(mesh, 31415)
!   call check_mesh_neighbours(mesh, 314)
!   call check_mesh_neighbours(mesh, 9082)
 

  call nullify_mesh_object()
  call fin()


  contains

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

    call assert_eq(vol_global, 31.006280191406404d0, "Mesh volume not matching")
    call assert_eq(nneg_vol_global, 0, "Cells with negative volume detected")

  end subroutine

end program test_read_mesh
