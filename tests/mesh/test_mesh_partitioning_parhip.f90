!> Test that partitions a very simple graph using ParHIP
!
!  Sample graph - adapted from ParMETIS manual to use 1-indexing
! 
!     1 -- 2 -- 3 -- 4 -- 5
!     |    |    |    |    |
!     6 -- 7 -- 8 -- 9 --10
!     |    |    |    |    |
!    11 --12 --13 --14 --15
!

program test_mesh_partitioning_parhip

  use MPI

  use testing_lib
  use partitioning, only: partition_kway
  use kinds, only: ccs_int, ccs_long
  use types, only: topology

  implicit none

  type(topology) :: topo

  call init()
  call initialise_test()

  ! Partition
  allocate(topo%global_partition(topo%global_num_cells))
  call partition_kway(par_env, topo)
  
  if(par_env%proc_id == 0) then
     print *, topo%global_partition
  end if

  ! Check mesh after partitioning
  call check_global_cell_count()

  call clean_up()
  call fin()

contains

  subroutine initialise_test

    topo%global_num_cells = 15
  
    select type(par_env)
    type is (parallel_environment_mpi)
  
      if(par_env%num_procs == 3) then
  
        allocate(topo%local_partition(5))
        allocate(topo%xadj(6))
        allocate(topo%vwgt(5)) 
        allocate(topo%vtxdist(4))
  
        if(par_env%proc_id == 0) then
          allocate(topo%adjncy(13))
          allocate(topo%adjwgt(13))
          topo%xadj = (/ 1, 3, 6, 9, 12, 14 /)
          topo%adjncy = (/ 2, 6, 1, 3, 7, 2, 4, 8, 3, 5, 9, 4, 10 /)
        else if (par_env%proc_id == 1) then
          allocate(topo%adjncy(18))
          allocate(topo%adjwgt(18))
          topo%xadj = (/ 1, 4, 8, 12, 16, 19 /)      
          topo%adjncy = (/ 1, 7, 11, 2, 6, 8, 12, 3, 7, 9, 13, 4, 8, 10, 14, 5, 9, 15 /)
        else 
          allocate(topo%adjncy(13))
          allocate(topo%adjwgt(13))
          topo%xadj = (/ 1, 3, 6, 9, 12, 14 /)
          topo%adjncy = (/ 6, 12, 7, 11, 13, 8, 12, 14, 9, 13, 15, 10, 14 /)
        end if
  
        topo%vtxdist = (/ 1, 6, 11, 16 /)
        topo%adjwgt = 1
        topo%vwgt = 1
  
      else
        write(message, *) "Test must be run on 3 MPI ranks"
        call stop_test(message)
      end if 
   
      class default
        write(message, *) "ERROR: Unknown parallel environment!"
        call stop_test(message)
    end select

  end subroutine

  subroutine check_global_cell_count

    integer :: nnew
    integer :: ntotal
    integer :: i

    nnew = 0
    do i = 1, topo%global_num_cells
       if (topo%global_partition(i) == par_env%proc_id) then
          nnew = nnew + 1
       end if
    end do
    call MPI_Allreduce(nnew, ntotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (ntotal /= topo%global_num_cells) then
       write(message, *) "ERROR: Total cell count after partitioning = ", ntotal, " expected ", topo%global_num_cells
       call stop_test(message)
    end if

  end subroutine

  subroutine clean_up

    if(allocated(topo%xadj)) then
      deallocate(topo%xadj)
    end if
  
    if(allocated(topo%adjncy)) then
      deallocate(topo%adjncy)
    end if
  
    if(allocated(topo%adjwgt)) then
      deallocate(topo%adjwgt)
    end if
  
    if(allocated(topo%vwgt)) then
      deallocate(topo%vwgt)
    end if
  
    if(allocated(topo%vtxdist)) then
      deallocate(topo%vtxdist)
    end if

  end subroutine

end program
