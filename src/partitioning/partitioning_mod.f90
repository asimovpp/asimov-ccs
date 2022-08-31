!> Module file partitioning.mod
!>
!> Module defining the partitioning interface for ASiMoV-CCS

module partitioning

  use kinds, only: ccs_long
  use types, only: topology
  use parallel_types, only: parallel_environment

  implicit none

  private
  public :: partition_kway
  public :: compute_partitioner_input
  ! public :: read_topology
  public :: compute_connectivity

  interface 

    !v Partition the mesh
    module subroutine partition_kway(par_env, topo)
      class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
      type(topology), target, intent(inout) :: topo                           !< The topology for which to compute the parition
    end subroutine partition_kway

    !v Compute the input arrays for the partitioner
    module subroutine compute_partitioner_input(par_env, topo)
      class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
      type(topology), target, intent(inout) :: topo                           !< The topology for which to compute the parition
    end subroutine compute_partitioner_input

    ! !v Read the topology data from an input (HDF5) file
    ! module subroutine read_topology(par_env, case_name, topo)
    !   class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    !   character(len=:), allocatable :: case_name
    !   type(topology), target, intent(inout) :: topo                           !< The mesh topology for which to compute the partition
    ! end subroutine read_topology

    ! !v Read the geometry data from an input (HDF5) file
    ! module subroutine read_geometry(par_env, case_name, geo)
    !   class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    !   character(len=:), allocatable :: case_name
    !   type(geometry), target, intent(inout) :: geo                            !< The mesh geometry
    ! end subroutine read_geometry

    module subroutine compute_connectivity(par_env, topo)
      class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
      type(topology), target, intent(inout) :: topo                           !< The topology for which to compute the parition
    end subroutine compute_connectivity
 
    end interface

end module partitioning