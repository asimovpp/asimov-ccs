submodule(partitioning) partitioning_parhip
#include "ccs_macros.inc"

  use kinds, only: ccs_int, ccs_long
  use types, only: topology, graph_connectivity
  use utils, only: str, debug_print
  use parallel_types_mpi, only: parallel_environment_mpi
  use meshing, only: get_global_num_cells
  use parallel, only: is_root, is_valid, create_shared_array, destroy_shared_array, sync
  use logging, only: log_unit_out
 
  implicit none

  interface
    subroutine partition_parhipkway(vtxdist, xadj, adjncy, vwgt, adjwgt, &
                                    num_procs, imbalance, suppress, &
                                    seed, mode, edgecuts, local_partition, comm) bind(c)
      use iso_c_binding

      integer(c_long), dimension(*) :: vtxdist
      integer(c_long), dimension(*) :: xadj
      integer(c_long), dimension(*) :: adjncy
      integer(c_long), dimension(*) :: vwgt
      integer(c_long), dimension(*) :: adjwgt
      integer(c_int) :: num_procs
      real(c_double) :: imbalance
      integer(c_int) :: suppress
      integer(c_int) :: seed
      integer(c_int) :: mode
      integer(c_int) :: edgecuts
      integer(c_long), dimension(*) :: local_partition
      integer(c_int) :: comm
    end subroutine
  end interface

contains

  !v Partition the mesh
  !
  ! Use ParHIP library to compute a k-way vertex separator given a k-way partition of the graph.
  ! The graph can be weighted or unweighted.
  !
  ! High-level interface operating on the mesh object.
  module subroutine partition_kway_mesh(par_env, mesh)

    class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The global parallel environment
    type(ccs_mesh), target, intent(inout) :: mesh                              !< The mesh for which to compute the parition

    call partition_kway_topo(par_env, mesh%topo)
    
  end subroutine partition_kway_mesh

  !v Partition the mesh
  !
  ! Use ParHIP library to compute a k-way vertex separator given a k-way partition of the graph.
  ! The graph can be weighted or unweighted.
  !
  ! High-level interface operating on the topology object.
  module subroutine partition_kway_topo(par_env, topo)

    class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The global parallel environment
    type(topology), target, intent(inout) :: topo                              !< The mesh topology for which to compute the parition

    call partition_kway_graph_conn(par_env, topo%graph_conn)
    
  end subroutine partition_kway_topo

  !v Partition the mesh
  !
  ! Use ParHIP library to compute a k-way vertex separator given a k-way partition of the graph.
  ! The graph can be weighted or unweighted.
  !
  ! Performs the partitioning on the graph connectivity object.
  module subroutine partition_kway_graph_conn(par_env, graph_conn)

    use mpi
    use iso_c_binding

    class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The global parallel environment
    type(graph_connectivity), target, intent(inout) :: graph_conn              !< The graph connectivity for which to compute the parition

    ! Local variables
    integer(ccs_int) :: local_part_size
    integer(ccs_int) :: irank

    integer(ccs_long) :: global_num_cells
    
    integer(c_long), dimension(:), allocatable :: vtxdist
    integer(c_long), dimension(:), allocatable :: xadj
    integer(c_long), dimension(:), allocatable :: adjncy
    integer(c_long), dimension(:), allocatable :: vwgt
    integer(c_long), dimension(:), allocatable :: adjwgt
    integer(c_int) :: num_procs
    real(c_double) :: imbalance
    integer(c_int) :: suppress
    integer(c_int) :: seed
    integer(c_int) :: mode
    integer(c_int) :: edgecuts
    integer(c_long), dimension(:), allocatable :: local_partition
    integer(c_int) :: comm

    ! Values hardcoded for now
    imbalance = 0.03  ! Desired balance - 0.03 = 3%
    seed = 2022       ! "Random" seed
    mode = 4          ! FASTSOCIAL
    suppress = 0      ! Do not suppress the output
    edgecuts = -1     ! XXX: silence unused variable warning

    irank = par_env%proc_id ! Current rank
    num_procs = par_env%num_procs

    ! vtxdist should contain the (initial) global cell partition - the final element is global cell
    ! count + 1
    global_num_cells = graph_conn%vtxdist(num_procs + 1) - 1
    
    ! ParHIP needs 0-indexing - shift array contents by -1
    allocate(vtxdist, mold=graph_conn%vtxdist)
    vtxdist = graph_conn%vtxdist - 1
    xadj = graph_conn%xadj - 1
    adjncy = graph_conn%adjncy - 1

    adjwgt = graph_conn%adjwgt
    vwgt = graph_conn%vwgt

    ! Set weights to 1
    adjwgt = 1_ccs_long
    vwgt = 1_ccs_long

    ! Number of elements in local partition array
    ! Needed for gathering loca partitions into global partition array
    local_part_size = size(graph_conn%local_partition)

    allocate (local_partition(local_part_size))

    ! Partitioning an unweighted graph
    select type (par_env)
    type is (parallel_environment_mpi)

      if (is_root(par_env)) then
        write(log_unit_out,*) "Partitioning with ParHIP"
      end if

      comm = par_env%comm

      call partition_parhipkway(vtxdist, xadj, adjncy, vwgt, adjwgt, &
                                num_procs, imbalance, suppress, &
                                seed, mode, edgecuts, local_partition, comm)

      graph_conn%local_partition(:) = local_partition(:)

    class default
      write(log_unit_out,*) "ERROR: Unknown parallel environment!"
    end select

    call dprint("Number of edgecuts: " // str(edgecuts))

  end subroutine partition_kway_graph_conn

  !v Compute the input arrays for the partitioner
  !
  ! Using the topology object, compute the input arrays for the ParHIP partitioner
  ! Input arrays for the partitioner are: vtxdist, xadj and adjncy
  module subroutine compute_partitioner_input(par_env, shared_env, mesh)

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The parallel environment
    type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the parition

    call compute_partitioner_input_generic(par_env, shared_env, mesh)

  end subroutine compute_partitioner_input

end submodule
