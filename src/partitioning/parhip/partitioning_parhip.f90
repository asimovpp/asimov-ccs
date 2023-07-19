submodule(partitioning) partitioning_parhip
#include "ccs_macros.inc"

  use kinds, only: ccs_int, ccs_real, ccs_long
  use utils, only: str, debug_print
  use parallel_types_mpi, only: parallel_environment_mpi
  use meshing, only: set_local_num_cells, get_local_num_cells, get_global_num_cells, &
                     get_max_faces

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
  module subroutine partition_kway(par_env, mesh)

    use mpi
    use iso_c_binding

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the parition

    ! Local variables
    integer(ccs_long), dimension(:), pointer :: tmp_partition
    integer(ccs_int) :: tmp_partition_window
    integer(ccs_int) :: local_part_size
    integer(ccs_int) :: irank
    integer(ccs_int) :: ierr
    integer(ccs_int) :: i

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

    integer(ccs_int) :: global_num_cells

    ! Values hardcoded for now
    imbalance = 0.03  ! Desired balance - 0.03 = 3%
    seed = 2022       ! "Random" seed
    mode = 4          ! FASTSOCIAL
    suppress = 0      ! Do not suppress the output
    edgecuts = -1     ! XXX: silence unused variable warning

    call get_global_num_cells(mesh, global_num_cells)

    irank = par_env%proc_id ! Current rank
    num_procs = par_env%num_procs

    ! ParHIP needs 0-indexing - shift array contents by -1
    vtxdist = mesh%topo%vtxdist - 1
    xadj = mesh%topo%xadj - 1
    adjncy = mesh%topo%adjncy - 1

    adjwgt = mesh%topo%adjwgt
    vwgt = mesh%topo%vwgt

    ! Set weights to 1
    adjwgt = 1
    vwgt = 1

    ! Number of elements in local partition array
    ! Needed for gathering loca partitions into global partition array
    local_part_size = size(mesh%topo%local_partition)

    allocate (local_partition(local_part_size))

    ! Partitioning an unweighted graph
    select type (par_env)
    type is (parallel_environment_mpi)

      comm = par_env%comm

      call partition_parhipkway(vtxdist, xadj, adjncy, vwgt, adjwgt, &
                                num_procs, imbalance, suppress, &
                                seed, mode, edgecuts, local_partition, comm)

      mesh%topo%local_partition(:) = local_partition(:)

      if (isroot(shared_env)) then
        call create_shared_array(shared_env, global_num_cells, tmp_partition, tmp_partition_window)
        tmp_partition(:) = 0
      else
        call create_shared_array(shared_env, 0, tmp_partition, tmp_partition_window)
      end if

      do i = 1, local_part_size
        tmp_partition(i + vtxdist(irank + 1)) = mesh%topo%local_partition(i)
      end do

      if (isroot(shared_env)) then
        call MPI_AllReduce(tmp_partition, mesh%topo%global_partition, global_num_cells, &
                          MPI_LONG, MPI_SUM, partition_env%comm, ierr)
      end if

    class default
      print *, "ERROR: Unknown parallel environment!"
    end select

    call dprint("Number of edgecuts: " // str(edgecuts))

    deallocate (tmp_partition)

  end subroutine partition_kway

  !v Compute the input arrays for the partitioner
  !
  ! Using the topology object, compute the input arrays for the ParHIP partitioner
  ! Input arrays for the partitioner are: vtxdist, xadj and adjncy
  module subroutine compute_partitioner_input(par_env, mesh)

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the parition

    call compute_partitioner_input_generic(par_env, mesh)

  end subroutine compute_partitioner_input

end submodule
