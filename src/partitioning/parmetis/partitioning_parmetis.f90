submodule(partitioning) partitioning_parmetis
#include "ccs_macros.inc"

  use kinds, only: ccs_int, ccs_long
  use types, only: topology, graph_connectivity
  use utils, only: str, debug_print
  use parallel_types_mpi, only: parallel_environment_mpi
  use meshing, only: set_local_num_cells, get_local_num_cells, get_global_num_cells
  use parallel, only: is_root, is_valid, create_shared_array, destroy_shared_array, sync
  use logging, only: log_unit_out

  implicit none

  interface
    subroutine partition_parmetiskway(vtxdist, xadj, adjncy, vwgt, adjwgt, &
                                      wgtflag, numflag, ncon, num_procs, &
                                      tpwgts, ubvec, options, &
                                      edgecuts, local_partition, comm) bind(c)
      use iso_c_binding

      integer(c_int32_t), dimension(*) :: vtxdist
      integer(c_int32_t), dimension(*) :: xadj
      integer(c_int32_t), dimension(*) :: adjncy
      integer(c_int32_t), dimension(*) :: vwgt
      integer(c_int32_t), dimension(*) :: adjwgt
      integer(c_int32_t) :: wgtflag ! Set to 0 for "no weights"
      integer(c_int32_t) :: numflag ! Numbering scheme - 1 means Fortran style
      integer(c_int32_t) :: ncon
      integer(c_int32_t) :: num_procs
      real(c_float), dimension(*) :: tpwgts
      real(c_float), dimension(*) :: ubvec
      integer(c_int32_t), dimension(*) :: options
      integer(c_int32_t) :: edgecuts
      integer(c_int32_t), dimension(*) :: local_partition
      integer(c_int) :: comm
    end subroutine
  end interface

contains

  !v Partition the mesh
  !
  ! Use Parmetis library to compute a k-way vertex separator given a k-way partition of the graph.
  ! The graph can be weighted or unweighted.
  !
  ! High-level interface operating on the mesh object.
  module subroutine partition_kway_mesh(par_env, mesh)

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the parition

    call partition_kway_topo(par_env, mesh%topo)
    
  end subroutine partition_kway_mesh

  !v Partition the mesh
  !
  ! Use Parmetis library to compute a k-way vertex separator given a k-way partition of the graph.
  ! The graph can be weighted or unweighted.
  !
  ! High-level interface operating on the topology object.
  module subroutine partition_kway_topo(par_env, topo)

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(topology), target, intent(inout) :: topo                           !< The mesh topology for which to compute the parition

    call partition_kway_graph_conn(par_env, topo%graph_conn)

  end subroutine partition_kway_topo

  !v Partition the mesh
  !
  ! Use Parmetis library to compute a k-way vertex separator given a k-way partition of the graph.
  ! The graph can be weighted or unweighted.
  !
  ! Performs the partitioning on the graph connectivity object.
  module subroutine partition_kway_graph_conn(par_env, graph_conn)

    use mpi
    use iso_c_binding
    use iso_fortran_env

    class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The parallel environment
    type(graph_connectivity), target, intent(inout) :: graph_conn              !< The graph connectivity for which to compute the parition

    ! Local variables
    integer(ccs_int) :: local_part_size
    integer(ccs_int) :: irank

    integer(ccs_long) :: global_num_cells

    integer(c_int32_t), dimension(:), allocatable :: vtxdist
    integer(c_int32_t), dimension(:), allocatable :: xadj
    integer(c_int32_t), dimension(:), allocatable :: adjncy
    integer(c_int32_t), dimension(:), allocatable :: vwgt
    integer(c_int32_t), dimension(:), allocatable :: adjwgt
    integer(c_int32_t) :: wgtflag ! Set to 0 for "no weights"
    integer(c_int32_t) :: numflag ! Numbering scheme - 1 means Fortran style
    integer(c_int32_t) :: ncon
    integer(c_int32_t) :: num_procs
    real(c_float), dimension(:), allocatable :: tpwgts
    real(c_float), dimension(:), allocatable :: ubvec
    integer(c_int32_t), dimension(:), allocatable :: options
    integer(c_int32_t) :: edgecuts
    integer(c_int32_t), dimension(:), allocatable :: local_partition
    integer(c_int) :: comm

    ! Values mostly hardcoded for now
    wgtflag = 0 ! No weights
    numflag = 0 ! Use C-style indexing for now
    ncon = 1
    num_procs = par_env%num_procs
    edgecuts = -1

    ! vtxdist should contain the (initial) global cell partition - the final element is global cell
    ! count + 1
    global_num_cells = graph_conn%vtxdist(num_procs + 1) - 1

    allocate (ubvec(ncon))
    allocate (tpwgts(ncon * num_procs))
    allocate (options(0:2))

    options(0) = 1 ! 0 = default values, 1 = values specified in (1) and (2)
    options(1) = 2 ! Output verbosity - 1 gives timing information
    options(2) = 2023 ! Random number seed

    ubvec(:) = 1.03 ! Imbalance tolerance for each vertex weight, 1.05 is recommended value
    tpwgts(:) = 1.0 / real(num_procs, c_float) ! Should be set to 1 / num_procs

    irank = par_env%proc_id ! Current rank

    vtxdist = int(graph_conn%vtxdist, int32) - 1
    xadj = int(graph_conn%xadj, int32) - 1
    adjncy = int(graph_conn%adjncy, int32) - 1

    adjwgt = int(graph_conn%adjwgt, int32)
    vwgt = int(graph_conn%vwgt, int32)

    ! Set weights to 1
    adjwgt = 1
    vwgt = 1

    ! Number of elements in local partition array
    ! Needed for gathering loca partitions into global partition array
    local_part_size = size(graph_conn%local_partition)

    allocate (local_partition(local_part_size))

    ! Partitioning an unweighted graph
    select type (par_env)
    type is (parallel_environment_mpi)

      if (is_root(par_env)) then
        write(log_unit_out,*) "Partitioning with ParMETIS"
      end if

      comm = par_env%comm

      call partition_parmetiskway(vtxdist, xadj, adjncy, vwgt, adjwgt, &
                                  wgtflag, numflag, ncon, num_procs, &
                                  tpwgts, ubvec, options, &
                                  edgecuts, local_partition, comm)

      graph_conn%local_partition(:) = int(local_partition(:), int64)

    class default
      write(log_unit_out,*) "ERROR: Unknown parallel environment! "
    end select

    call dprint("Number of edgecuts: " // str(int(edgecuts)))

  end subroutine partition_kway_graph_conn

  !v Compute the input arrays for the partitioner
  !
  ! Using the topology object, compute the input arrays for the Parmetis partitioner
  ! Input arrays for the partitioner are: vtxdist, xadj and adjncy
  module subroutine compute_partitioner_input(par_env, shared_env, mesh)

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The parallel environment
    type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the parition

    call compute_partitioner_input_generic(par_env, shared_env, mesh)

  end subroutine compute_partitioner_input

end submodule
