submodule(partitioning) partitioning_parmetis
#include "ccs_macros.inc"

  use kinds, only: ccs_int, ccs_real, ccs_long
  use types, only: topology, graph_connectivity
  use utils, only: str, debug_print
  use parallel_types_mpi, only: parallel_environment_mpi
  use meshing, only: set_local_num_cells, get_local_num_cells, get_global_num_cells
  use parallel, only: is_root, is_valid, create_shared_array, destroy_shared_array, sync

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
  module subroutine partition_kway(par_env, shared_env, roots_env, mesh)

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: roots_env !< The parallel environment
    type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the parition

    call partition_kway_topo(par_env, shared_env, roots_env, mesh%topo)
    
  end subroutine partition_kway

  !v Partition the mesh
  !
  ! Use Parmetis library to compute a k-way vertex separator given a k-way partition of the graph.
  ! The graph can be weighted or unweighted.
  !
  ! High-level interface operating on the topology object.
  subroutine partition_kway_topo(par_env, shared_env, roots_env, topo)

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: roots_env !< The parallel environment
    type(topology), target, intent(inout) :: topo                           !< The mesh topology for which to compute the parition

    integer(ccs_int) :: global_num_cells

    call get_global_num_cells(topo, global_num_cells)

    call partition_kway_graphconn(par_env, shared_env, roots_env, global_num_cells, topo%graph_conn)

  end subroutine partition_kway_topo

  !v Partition the mesh
  !
  ! Use Parmetis library to compute a k-way vertex separator given a k-way partition of the graph.
  ! The graph can be weighted or unweighted.
  !
  ! Performs the partitioning on the graph connectivity object.
  subroutine partition_kway_graphconn(par_env, shared_env, roots_env, global_num_cells, graph_conn)

    use mpi
    use iso_c_binding
    use iso_fortran_env

    class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: roots_env  !< The parallel environment
    integer(ccs_int), intent(in) :: global_num_cells                           !< The global cell count
    type(graph_connectivity), target, intent(inout) :: graph_conn              !< The graph connectivity for which to compute the parition

    ! Local variables
    integer(ccs_long), dimension(:), pointer :: tmp_partition
    integer(ccs_int) :: tmp_partition_window
    integer(ccs_int) :: local_part_size
    integer(ccs_int) :: irank
    integer(ccs_int) :: ierr
    integer(ccs_int) :: i

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
    ncon = 3
    num_procs = par_env%num_procs
    edgecuts = -1

    allocate (ubvec(ncon))
    allocate (tpwgts(ncon * num_procs))
    allocate (options(0:2))

    options(0) = 0 ! 0 = default values, 1 = values specified in (1) and (2)
    options(1) = 1 ! Output verbosity - 1 gives timing information
    options(2) = 2023 ! Random number seed

    ubvec(:) = 1.05 ! Imbalance tolerance for each vertex weight, 1.05 is recommended value
    tpwgts(:) = 1.0 / real(num_procs, c_float) ! Sum of tpwgts(:) should be 1. Check this is correct

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

      comm = par_env%comm

      call partition_parmetiskway(vtxdist, xadj, adjncy, vwgt, adjwgt, &
                                  wgtflag, numflag, ncon, num_procs, &
                                  tpwgts, ubvec, options, &
                                  edgecuts, local_partition, comm)

      graph_conn%local_partition(:) = int(local_partition(:), int64)

      call create_shared_array(shared_env, global_num_cells, tmp_partition, tmp_partition_window)

      if (is_root(shared_env)) then
        tmp_partition(:) = 0
      end if
      call sync(shared_env)

      do i = 1, local_part_size
        tmp_partition(i + vtxdist(irank + 1)) = graph_conn%local_partition(i)
      end do
      call sync(shared_env)

      if (is_valid(roots_env)) then
        select type (roots_env)
        type is (parallel_environment_mpi)
          call MPI_AllReduce(tmp_partition, graph_conn%global_partition, global_num_cells, &
                            MPI_INTEGER8, MPI_SUM, roots_env%comm, ierr)
        class default
          print *, "ERROR: Unknown parallel environment!"
        end select
      end if

    class default
      print *, "ERROR: Unknown parallel environment! "
    end select

    call dprint("Number of edgecuts: " // str(int(edgecuts)))

    call sync(shared_env)

    call destroy_shared_array(shared_env, tmp_partition, tmp_partition_window)

  end subroutine partition_kway_graphconn

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
