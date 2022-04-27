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

  use testing_lib
  use partitioning, only: partition_kway
  use kinds, only: ccs_int, ccs_long
  use types, only: topology

  type(topology) :: topo
  integer(ccs_long), dimension(:), allocatable :: partition_vector

  implicit none

  call init()

  select type(par_env)
  type is (parallel_environment_mpi)

    if(par_env%num_procs == 3) then

      if(par_env%proc_id == 0) then
        allocate(topo%xadj(6))
        allocate(topo%adjncy(13))
        allocate(topo%vtxdist(4))
        topo%xadj = (/ 1, 3, 6, 9, 12, 14 /)
        topo%adjncy = (/ 2, 6, 1, 3, 7, 2, 4, 8, 3, 5, 9, 4, 10 /)
        topo%vtxdist = (/ 1, 6, 11, 16 /)
      else if (par_env%proc_id == 1) then
        allocate(topo%xadj(6))
        allocate(topo%adjncy(15))
        allocate(topo%vtxdist(4))
        topo%xadj = (/ 1, 4, 8, 12, 16, 19 /)      
        topo%adjncy = (/ 1, 7, 11, 2, 6, 8, 12, 3, 7, 9, 13, 4, 8, 10, 14, 5, 9, 15 /)
        topo%vtxdist = (/ 1, 6, 11, 16 /)
      else 
        allocate(topo%xadj(6))
        allocate(topo%adjncy(13))
        allocate(topo%vtxdist(4))
        topo%xadj = (/ 1, 3, 6, 9, 12, 14 /)
        topo%adjncy = (/ 6, 12, 7, 11, 13, 8, 12, 14, 9, 13, 15, 10, 14 /)
        topo%vtxdist = (/ 1, 6, 11, 16 /)
      end if
    
    else
      print*,"Test must be run on 3 MPI ranks"
    end if 

    ! call partition_kway(par_env, topo, partition_vector)

  class default
    write(message, *) "ERROR: Unknown parallel environment!"
    call stop_test(message)
  end select

  if(allocated(topo%xadj)) then
    deallocate(topo%xadj)
  end if

  if(allocated(topo%adjncy)) then
    deallocate(topo%adjncy)
  end if

  if(allocated(topo%vtxdist)) then
    deallocate(topo%vtxdist)
  end if

  call fin()

end program