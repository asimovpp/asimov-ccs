submodule (partitioning) partitioning_parhip

  use kinds, only: ccs_int, ccs_real

  implicit none

contains

  !v Partition the mesh
  !
  ! Use ParHIP library to compute a k-way vertex separator given a k-way partition of the graph.
  ! The graph can be weighted or unweighted.
  module subroutine partition_kway(par_env, topo, partition_vector)

    use iso_c_binding

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(topology), target, intent(inout) :: topo                           !< The topology for which to compute the parition
    integer(ccs_long), allocatable, intent(out) :: partition_vector         !< The vector with the partition 

    ! Local variables
    real(ccs_real) :: imbalance
    integer(ccs_int) :: seed
    integer(ccs_int) :: mode
    integer(ccs_int) :: suppress
    type(c_ptr) :: vwgt     ! NULL pointer to be used if no weights are attached to vertices
    type(c_ptr) :: adjwgt   ! NULL pointer to be used if no weights are attached to edges

    ! Values hardcoded for now
    imbalance = 0.03  ! Desired balance - 0.03 = 3% 
    seed = 2022       ! "Random" seed
    mode = 3          ! FASTSOCIAL
    suppress = 0      ! Do not suppress the output

    ! ParHIP needs 0-indexing
    topo%vtxdist = topo%vtxdist -1
    topo%xadj = topo%xadj - 1
    topo%adjncy = topo%adjncy - 1

    if(allocated(topo%vwgt) .eqv. .false. .and. allocated(topo%adjwgt) .eqv. .false.) then

      vwgt = c_null_ptr
      adjwgt = c_null_ptr
   
      ! Partitioning an unweighted graph
      ! call partition_parhipkway(topo%vtxdist, topo%xadj, topo%adjncy, vwgt, adjwgt, & 
                                  ! par_env%num_procs, imbalance, suppress, seed, &
                                  ! mode, edgecut, partition_vector, par_env%comm)
      
    else 

      ! Partitioning a graph with weights on vertices and edges
      ! call partition_parhipkway(topo%vtxdist, topo%xadj, topo%adjncy, topo%vwgt, topo%adjwgt, & 
                                  ! par_env%num_procs, imbalance, suppress, seed, &
                                  ! mode, edgecut, partition_vector, par_env%comm)

    end if

    ! Return to 1-indexing. For safety only right now, later on these arrays
    ! will be deallocated as they are not used beyond this point
    topo%vtxdist = topo%vtxdist -1
    topo%xadj = topo%xadj - 1
    topo%adjncy = topo%adjncy - 1

    ! Make sure partition vector starts from 1
    partition_vector = partition_vector + 1

  end subroutine partition_kway

end submodule