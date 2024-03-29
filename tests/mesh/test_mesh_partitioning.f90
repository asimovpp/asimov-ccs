!v Test that partitions a very simple graph
!
!  Sample graph - adapted from ParMETIS manual to use 1-indexing
!
!  1 -- 2 -- 3 -- 4 -- 5
!  |    |    |    |    |
!  6 -- 7 -- 8 -- 9 --10
!  |    |    |    |    |
!  11 --12 --13 --14 --15

program test_mesh_partitioning

  use MPI

  use testing_lib
  use partitioning, only: partition_kway
  use kinds, only: ccs_int, ccs_long
  use types, only: ccs_mesh, topology
  use meshing, only: get_global_num_cells, set_global_num_cells
  use meshing, only: set_mesh_object, nullify_mesh_object

  implicit none

  integer(ccs_int) :: global_num_cells

  call init()

  call set_mesh_object(mesh)
  call initialise_test()

  ! Partition
  call get_global_num_cells(global_num_cells)
  call create_shared_array(shared_env, global_num_cells, mesh%topo%graph_conn%global_partition, mesh%topo%graph_conn%global_partition_window)
  call partition_kway(par_env, shared_env, roots_env, mesh)

  if (par_env%proc_id == 0) then
    print *, mesh%topo%graph_conn%global_partition
  end if

  ! Check mesh after partitioning
  call check_global_cell_count()

  call nullify_mesh_object()
  call clean_up()
  call fin()

contains

  subroutine initialise_test()

    call set_global_num_cells(15)

    select type (par_env)
    type is (parallel_environment_mpi)

      if (par_env%num_procs == 3) then

        allocate (mesh%topo%graph_conn%local_partition(5))
        allocate (mesh%topo%graph_conn%xadj(6))
        allocate (mesh%topo%graph_conn%vwgt(5))
        allocate (mesh%topo%graph_conn%vtxdist(4))

        if (par_env%proc_id == 0) then
          allocate (mesh%topo%graph_conn%adjncy(13))
          allocate (mesh%topo%graph_conn%adjwgt(13))
          mesh%topo%graph_conn%xadj = (/1, 3, 6, 9, 12, 14/)
          mesh%topo%graph_conn%adjncy = (/2, 6, 1, 3, 7, 2, 4, 8, 3, 5, 9, 4, 10/)
        else if (par_env%proc_id == 1) then
          allocate (mesh%topo%graph_conn%adjncy(18))
          allocate (mesh%topo%graph_conn%adjwgt(18))
          mesh%topo%graph_conn%xadj = (/1, 4, 8, 12, 16, 19/)
          mesh%topo%graph_conn%adjncy = (/1, 7, 11, 2, 6, 8, 12, 3, 7, 9, 13, 4, 8, 10, 14, 5, 9, 15/)
        else
          allocate (mesh%topo%graph_conn%adjncy(13))
          allocate (mesh%topo%graph_conn%adjwgt(13))
          mesh%topo%graph_conn%xadj = (/1, 3, 6, 9, 12, 14/)
          mesh%topo%graph_conn%adjncy = (/6, 12, 7, 11, 13, 8, 12, 14, 9, 13, 15, 10, 14/)
        end if

        mesh%topo%graph_conn%vtxdist = (/1, 6, 11, 16/)
        mesh%topo%graph_conn%adjwgt = 1
        mesh%topo%graph_conn%vwgt = 1

      else
        write (message, *) "Test must be run on 3 MPI ranks"
        call stop_test(message)
      end if

    class default
      write (message, *) "ERROR: Unknown parallel environment!"
      call stop_test(message)
    end select

  end subroutine

  subroutine check_global_cell_count

    integer :: nnew
    integer :: ntotal
    integer :: i

    integer(ccs_int) :: global_num_cells

    nnew = 0
    call get_global_num_cells(global_num_cells)
    do i = 1, global_num_cells
      if (mesh%topo%graph_conn%global_partition(i) == par_env%proc_id) then
        nnew = nnew + 1
      end if
    end do
    call MPI_Allreduce(nnew, ntotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (ntotal /= global_num_cells) then
      write (message, *) "ERROR: Total cell count after partitioning = ", ntotal, " expected ", global_num_cells
      call stop_test(message)
    end if

  end subroutine

  subroutine clean_up

    if (allocated(mesh%topo%graph_conn%xadj)) then
      deallocate (mesh%topo%graph_conn%xadj)
    end if

    if (allocated(mesh%topo%graph_conn%adjncy)) then
      deallocate (mesh%topo%graph_conn%adjncy)
    end if

    if (allocated(mesh%topo%graph_conn%adjwgt)) then
      deallocate (mesh%topo%graph_conn%adjwgt)
    end if

    if (allocated(mesh%topo%graph_conn%vwgt)) then
      deallocate (mesh%topo%graph_conn%vwgt)
    end if

    if (allocated(mesh%topo%graph_conn%vtxdist)) then
      deallocate (mesh%topo%graph_conn%vtxdist)
    end if

  end subroutine

end program
