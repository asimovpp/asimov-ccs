!> Module file partitioning.mod
!>
!> Module defining the partitioning interface for ASiMoV-CCS

module partitioning

  use kinds, only: ccs_long
  use types, only: ccs_mesh
  use parallel_types, only: parallel_environment

  implicit none

  private
  public :: partition_kway
  public :: compute_partitioner_input
  public :: compute_connectivity, compute_connectivity_get_local_cells

  interface

    !v Partition the mesh
    module subroutine partition_kway(par_env, mesh)
      class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
      type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the parition
    end subroutine partition_kway

    !v Compute the input arrays for the partitioner
    module subroutine compute_partitioner_input(par_env, mesh)
      class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
      type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the parition
    end subroutine compute_partitioner_input

    module subroutine compute_connectivity(par_env, mesh)
      class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
      type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the parition
    end subroutine compute_connectivity

    module subroutine compute_connectivity_get_local_cells(par_env, mesh)
      class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
      type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the parition
    end subroutine compute_connectivity_get_local_cells
    
  end interface

end module partitioning
