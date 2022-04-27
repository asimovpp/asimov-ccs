!> Test that partitions a very simple graph using ParHIP
!
!  Sample graph
! 
!   0 -- 1 -- 2 -- 3 -- 4
!   |    |    |    |    |
!   5 -- 6 -- 7 -- 8 -- 9
!   |    |    |    |    |
!  10 --11 --12 --13 --14
!

program test_mesh_partitioning_parhip

  use testing_lib
  use partitioning, only: partition_kway
  use kinds, only: ccs_int, ccs_idx
  use types, only: topology

  type(topology) :: topo
  integer(ccs_idx), dimension(:), allocatable :: partition_vector


  implicit none

  call init()

  select type(par_env)
  type is (parallel_environment_mpi)

    if(par_env%nprocs == 3) then

      if(par_env%proc_id == 0) then
        allocate(topo%xadj(6))
        allocate(topo%vtxdist(4))
        allocate(topo%adjncy(13))
        topo%xadj = (/ 0 2 5 8 11 13 /)
        topo%vtxdist = (/ 0 5 10 15 /)
        topo%adjncy = (/ 1 5 0 2 6 1 3 7 2 4 8 3 9 /)
      else if (par_env%proc_id == 1) then
        allocate(topo%xadj(6))
        allocate(topo%vtxdist(4))
        allocate(topo%adjncy(15))
        topo%xadj = (/ 0 3 7 11 15 18 /)      
        topo%vtxdist = (/ 0 5 10 15 /)
        topo%adjncy = (/ 0 6 10 1 5 7 11 2 6 8 12 3 7 9 13 4 8 14 /)
      else 
        allocate(topo%xadj(6))
        allocate(topo%vtxdist(4))
        allocate(topo%adjncy(13))
        topo%vtxdist = (/ 0 5 10 15 /)
        topo%xadj = (/ 0 2 5 8 11 13 /)
        topo%adjncy = (/ 5 11 6 10 12 7 11 13 8 12 14 9 13 /)
      end if
    
    else
      print*,"Test must be run on 3 MPI ranks"
    end if 

    !call partition_kway(par_env, topo, partition_vector)

  class default
    write(message, *) "ERROR: Unknown parallel environment!"
    call stop_test(message)
  end select

  if(allocated(adjncy)) then
    deallocate(adjncy)
  end if

  call fin()

end program