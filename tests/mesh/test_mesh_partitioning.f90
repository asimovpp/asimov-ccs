!v Test that partitions a very simple graph
!
!  Sample graph - adapted from ParMETIS manual to use 1-indexing
!
!  1 -- 2 -- 3 -- 4 -- 5
!  |    |    |    |    |
!  6 -- 7 -- 8 -- 9 --10
!  |    |    |    |    |
!  11 --12 --13 --14 --15

module m_test_mesh_partitioning

  use mpi

  use testing_lib, only: message, stop_test

  use kinds, only: ccs_int, ccs_long
  use types, only: graph_connectivity
  use parallel_types, only: parallel_environment

  implicit none

  private
  public :: initialise_test
  public :: check_global_cell_count
  public :: clean_up
  public :: get_global_num_cells

contains

  subroutine initialise_test(par_env, graph_conn)

    class(parallel_environment), intent(in) :: par_env
    type(graph_connectivity), intent(out) :: graph_conn

    if (par_env%num_procs == 3) then
       allocate (graph_conn%local_partition(5))
       allocate (graph_conn%xadj(6))
       allocate (graph_conn%vwgt(5))
       allocate (graph_conn%vtxdist(4))

       if (par_env%proc_id == 0) then
          allocate (graph_conn%adjncy(13))
          allocate (graph_conn%adjwgt(13))
          graph_conn%xadj = (/1, 3, 6, 9, 12, 14/)
          graph_conn%adjncy = (/2, 6, 1, 3, 7, 2, 4, 8, 3, 5, 9, 4, 10/)
       else if (par_env%proc_id == 1) then
          allocate (graph_conn%adjncy(18))
          allocate (graph_conn%adjwgt(18))
          graph_conn%xadj = (/1, 4, 8, 12, 16, 19/)
          graph_conn%adjncy = (/1, 7, 11, 2, 6, 8, 12, 3, 7, 9, 13, 4, 8, 10, 14, 5, 9, 15/)
       else
          allocate (graph_conn%adjncy(13))
          allocate (graph_conn%adjwgt(13))
          graph_conn%xadj = (/1, 3, 6, 9, 12, 14/)
          graph_conn%adjncy = (/6, 12, 7, 11, 13, 8, 12, 14, 9, 13, 15, 10, 14/)
       end if

       graph_conn%vtxdist = (/1, 6, 11, 16/)
       graph_conn%adjwgt = 1
       graph_conn%vwgt = 1
    else
       write (message, *) "Test must be run on 3 MPI ranks"
       call stop_test(message)
    end if

  end subroutine

  subroutine check_global_cell_count(par_env, graph_conn)

    class(parallel_environment), intent(in) :: par_env
    type(graph_connectivity), intent(in) :: graph_conn

    integer, dimension(:), allocatable :: nnew
    integer :: ntotal

    integer(ccs_long) :: local_num_cells
    integer(ccs_long) :: i

    integer :: ierr

    allocate(nnew(par_env%num_procs))
    nnew = 0
    call get_local_num_cells(graph_conn, local_num_cells)
    do i = 1, local_num_cells
      nnew(graph_conn%local_partition(i) + 1) = nnew(graph_conn%local_partition(i) + 1) + 1
    end do
    call MPI_Allreduce(MPI_IN_PLACE, nnew, par_env%num_proces, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
    ntotal = nnew(par_env%proc_id + 1)
    if (ntotal /= global_num_cells) then
      write (message, *) "ERROR: Total cell count after partitioning = ", ntotal, " expected ", global_num_cells
      call stop_test(message)
    end if

  end subroutine

  subroutine clean_up(graph_conn)

    type(graph_connectivity), intent(inout) :: graph_conn

    if (allocated(graph_conn%xadj)) then
      deallocate (graph_conn%xadj)
    end if

    if (allocated(graph_conn%adjncy)) then
      deallocate (graph_conn%adjncy)
    end if

    if (allocated(graph_conn%adjwgt)) then
      deallocate (graph_conn%adjwgt)
    end if

    if (allocated(graph_conn%vwgt)) then
      deallocate (graph_conn%vwgt)
    end if

    if (allocated(graph_conn%vtxdist)) then
      deallocate (graph_conn%vtxdist)
    end if

  end subroutine clean_up

  subroutine get_global_num_cells(graph_conn, global_num_cells)

    type(graph_connectivity), intent(in) :: graph_conn
    integer(ccs_long), intent(out) :: global_num_cells

    integer(ccs_int) :: np

    np = size(graph_conn%vtxdist)
    global_num_cells = graph_conn%vtxdist(np) - 1

  end subroutine get_global_num_cells
  
end module m_test_mesh_partitioning

program test_mesh_partitioning

  use MPI

  use testing_lib
  use partitioning, only: partition_kway
  use kinds, only: ccs_int, ccs_long
  use types, only: ccs_mesh, topology, graph_connectivity

  use m_test_mesh_partitioning

  implicit none

  integer(ccs_long) :: global_num_cells
  type(graph_connectivity) :: graph_conn

  call init()

  call initialise_test(par_env, graph_conn)

  ! Partition
  call get_global_num_cells(graph_conn, global_num_cells)
  call create_shared_array(shared_env, int(global_num_cells, ccs_int), graph_conn%global_partition, graph_conn%global_partition_window)
  call partition_kway(par_env, shared_env, roots_env, graph_conn)

  if (par_env%proc_id == 0) then
    print *, graph_conn%global_partition
  end if

  ! Check mesh after partitioning
  call check_global_cell_count(par_env, graph_conn)

  call clean_up(graph_conn)
  call fin()

end program test_mesh_partitioning
