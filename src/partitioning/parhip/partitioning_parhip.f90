submodule (partitioning) partitioning_parhip

  implicit none

contains

  !v Partition the mesh
  !
  ! Use ParHIP library to compute a k-way vertex separator given a k-way partition of the graph
  module subroutine partition_kway(par_env, topo, partition_vector)

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(topology), target, intent(in) :: topo                              !< The topology for which to compute the parition
    integer(ccs_idx), allocatable, intent(out) :: partition_vector          !< The vector with the partition 

    ! Local variables
    real(accs_real) :: imbalance
    integer(accs_int) :: seed
    integer(accs_int) :: mode
    logical :: suppress

    ! TODO 
    ! call partition_parhipkway(vtxdist, xadj, adjncy, vwgt, adjwgt, & 
                                ! par_env%num_procs, imbalance, suppress, seed, &
                                ! mode, edgecut, partition_vector, par_env%comm)
    
  end subroutine partition_kway

end submodule