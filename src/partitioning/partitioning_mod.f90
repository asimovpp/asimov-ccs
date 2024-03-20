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
  public :: cleanup_partitioner_data
  public :: compute_connectivity
  public :: compute_connectivity_get_local_cells
  public :: print_partition_quality
  
  interface

    !v Partition the mesh
    module subroutine partition_kway(par_env, shared_env, roots_env, mesh)
      class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The global parallel environment
      class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared parallel environment
      class(parallel_environment), allocatable, target, intent(in) :: roots_env  !< The roots of shared memory parallel environment
      type(ccs_mesh), target, intent(inout) :: mesh                              !< The mesh for which to compute the partition
    end subroutine partition_kway

    !v Compute the input arrays for the partitioner
    module subroutine compute_partitioner_input(par_env, shared_env, mesh)
      class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The global parallel environment
      class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared parallel environment
      type(ccs_mesh), target, intent(inout) :: mesh                              !< The mesh for which to compute the partition
    end subroutine compute_partitioner_input
   
  !v Deallocate partitioner data structures associated with the mesh
    module subroutine cleanup_partitioner_data(shared_env, mesh)
      class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared parallel environment
      type(ccs_mesh), target, intent(inout) :: mesh                              !< The mesh
    end subroutine cleanup_partitioner_data

    module subroutine compute_connectivity(par_env,shared_env, mesh)
      class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The global parallel environment
      class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared parallel environment
      type(ccs_mesh), target, intent(inout) :: mesh                              !< The mesh for which to compute the partition
    end subroutine compute_connectivity

    module subroutine compute_connectivity_get_local_cells(par_env, mesh)
      class(parallel_environment), allocatable, target, intent(in) :: par_env !< The global parallel environment
      type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the partition
    end subroutine compute_connectivity_get_local_cells

    !v Internal routine for computingd the input arrays for the partitioner
    module subroutine compute_partitioner_input_generic(par_env, shared_env, mesh)
      class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The global parallel environment
      class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared parallel environment
      type(ccs_mesh), target, intent(inout) :: mesh                              !< The mesh for which to compute the partition
    end subroutine compute_partitioner_input_generic

    !v Compute and report the partitioning quality.
    !
    !  The following metrics are implemented
    !  - The "surface to volume ratio" nhalo / nlocal (averaged)
    !  - The minimum departure from load balance min(nlocal) / avg(nlocal)
    !  - The maximum departure from load balance max(nlocal) / avg(nlocal)
    module subroutine print_partition_quality(par_env) 
      class(parallel_environment), intent(in) :: par_env
    end subroutine print_partition_quality
    
  end interface

end module partitioning
