!> Module file partitioning.mod
!>
!> Module defining the partitioning interface for ASiMoV-CCS

module partitioning

  use core, only: ccs_options
  use kinds, only: ccs_int, ccs_long
  use types, only: ccs_mesh, topology, graph_connectivity
  use parallel_types, only: parallel_environment

  implicit none

  private
  public :: partition_kway
  public :: compute_partitioner_input
  public :: cleanup_partitioner_data
  public :: compute_connectivity
  public :: compute_connectivity_get_local_cells
  public :: print_partition_quality
  public :: proccnt_to_vtxdist
  public :: partition_count
  public :: compute_global_indices_partition

  interface partition_kway
    module procedure :: partition_kway_mesh
    module procedure :: partition_kway_topo
    module procedure :: partition_kway_graph_conn
  end interface partition_kway

  interface

    !v Partition the mesh
    module subroutine partition_kway_mesh(par_env, mesh)
      class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The global parallel environment
      type(ccs_mesh), target, intent(inout) :: mesh                              !< The mesh for which to compute the partition
    end subroutine partition_kway_mesh
    module subroutine partition_kway_topo(par_env, topo)
      class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The global parallel environment
      type(topology), target, intent(inout) :: topo                              !< The topo for which to compute the partition
    end subroutine partition_kway_topo
    module subroutine partition_kway_graph_conn(par_env, graph_conn)
      class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The global parallel environment
      type(graph_connectivity), target, intent(inout) :: graph_conn              !< The graph_conn for which to compute the partition
    end subroutine partition_kway_graph_conn

    !v Compute the input arrays for the partitioner
    module subroutine compute_partitioner_input(par_env, shared_env, mesh)
      class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The global parallel environment
      class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared parallel environment
      type(ccs_mesh), target, intent(inout) :: mesh                              !< The mesh for which to compute the partition
    end subroutine compute_partitioner_input

    !v Deallocate partitioner data structures associated with the mesh
    module subroutine cleanup_partitioner_data(mesh)
      type(ccs_mesh), target, intent(inout) :: mesh                              !< The mesh
    end subroutine cleanup_partitioner_data

    module subroutine compute_connectivity(par_env, shared_env, mesh)
      class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The global parallel environment
      class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared parallel environment
      type(ccs_mesh), target, intent(inout) :: mesh                              !< The mesh for which to compute the partition
    end subroutine compute_connectivity

    module subroutine compute_connectivity_get_local_cells(par_env, topo)
      class(parallel_environment), allocatable, target, intent(in) :: par_env !< The global parallel environment
      type(topology), target, intent(inout) :: topo                           !< The mesh topology for which to compute the partition
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
    module subroutine print_partition_quality(par_env, run_options)
      class(parallel_environment), intent(in) :: par_env
      type(ccs_options), intent(in) :: run_options
    end subroutine print_partition_quality

    !!! INTERNAL SUBROUTINES - ONLY DECLARED FOR TESTING !!!
    pure module function partition_count(nproc, partition) result(proc_count)
      integer(ccs_int), intent(in) :: nproc
      integer(ccs_long), dimension(:), intent(in) :: partition
      integer(ccs_int), dimension(nproc) :: proc_count
    end function partition_count
    module function proccnt_to_vtxdist(proccnt) result(vtxdist)
      integer(ccs_int), dimension(:), intent(in) :: proccnt
      integer(ccs_long), dimension(:), allocatable :: vtxdist
    end function proccnt_to_vtxdist
    module function compute_global_indices_partition(partition, proc_ctr, vtxdist, global_idx_start) result(global_indices)
      integer(ccs_long), dimension(:), intent(in) :: partition
      integer(ccs_int), dimension(:), intent(in) :: proc_ctr
      integer(ccs_long), dimension(:), intent(in) :: vtxdist
      integer(ccs_long), intent(in) :: global_idx_start
      integer(ccs_long), dimension(:), allocatable :: global_indices
    end function compute_global_indices_partition
    !!! INTERNAL SUBROUTINES - ONLY DECLARED FOR TESTING !!!

  end interface

end module partitioning
