program test_mesh_accessors

  use mpi

  use testing_lib
  use meshing, only: get_local_num_cells
  use mesh_utils, only: build_mesh
  
  implicit none

  integer(ccs_int) :: n, nx, ny, nz
  real(ccs_real) :: l
  type(ccs_mesh) :: mesh
  
  call init()
  
  l = 1.0_ccs_real
  nx = 5; ny = 5; nz = 5
  n = nx * ny * nz
  mesh = build_mesh(par_env, nx, ny, nz, l)

  call check_local_num_cells()

  call fin()

contains

  subroutine check_local_num_cells

    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: global_num_cells
    integer(ccs_err) :: ierr
    integer(ccs_long) :: local_num_cells_long

    character(len=256) :: message
    
    call get_local_num_cells(mesh, local_num_cells)
    if (local_num_cells > n) then
       print *, "FAIL: local_num_cells exceeds global cell count!"
       stop 1
    end if

    ! Check global sum
    select type(par_env)
    type is (parallel_environment_mpi)
       call MPI_Allreduce(local_num_cells, global_num_cells, 1, MPI_INT, MPI_SUM, par_env%comm, ierr)
       if (global_num_cells /= n) then
          write(message, *) "FAIL: global cell count is wrong!"
          call stop_test(message)
       end if
    class default
       write(message, *) "FAIL: Unknonwn parallel environment!"
       call stop_test(message)
    end select

    ! Check using long data
    call get_local_num_cells(mesh, local_num_cells_long)
    local_num_cells = int(local_num_cells, ccs_int) ! Given mesh size this should not exceed a normal int
    if (local_num_cells > n) then
       print *, "FAIL: local_num_cells exceeds global cell count!"
       stop 1
    end if

    ! Check global sum
    select type(par_env)
    type is (parallel_environment_mpi)
       call MPI_Allreduce(local_num_cells, global_num_cells, 1, MPI_INT, MPI_SUM, par_env%comm, ierr)
       if (global_num_cells /= n) then
          write(message, *) "FAIL: global cell count is wrong!"
          call stop_test(message)
       end if
    class default
       write(message, *) "FAIL: Unknonwn parallel environment!"
       call stop_test(message)
    end select
    
  end subroutine check_local_num_cells
  
end program test_mesh_accessors
