!> Module file partitioning.mod
!>
!> Module defining the partitioning interface for ASiMoV-CCS

module partitioning

  use kinds, only: ccs_idx
  use types, only: topology
  use parallel, only: parallel_environment

  implicit none

  private
  public :: partition_kway

  !v Partition the mesh
  module subroutine partition_kway(par_env, topo, partition_vector)
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(topology), target, intent(in) :: topo                              !< The topology for which to compute the parition
    integer(ccs_idx), allocatable, intent(out) :: partition_vector          !< The vector with the partition 
  end subroutine partition_kway

end module partitioning