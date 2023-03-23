module parhip_interfaces
#include "ccs_macros.inc"

    use kinds, only: ccs_int, ccs_real, ccs_long
    use utils, only: str, debug_print
    use parallel_types_mpi, only: parallel_environment_mpi
    use meshing, only: set_local_num_cells, get_local_num_cells
    use kinds, only: ccs_long
    use types, only: ccs_mesh
    use parallel_types, only: parallel_environment
  
  
    use iso_c_binding

  implicit none

  public 

  interface 
    subroutine partition_parhipkway(vtxdist, xadj, adjncy, vwgt, adjwgt, &
                                  num_procs, imbalance, suppress, &
                                  seed, mode, edgecuts, local_partition, comm) bind(c)
      use iso_c_binding

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
    end subroutine 
  end interface

  contains

  subroutine partition_kway_foo(par_env, mesh)

    use mpi
    use iso_c_binding

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the parition

    ! Local variables
    integer(ccs_long), dimension(:), allocatable :: tmp_partition
    ! real(ccs_real) :: imbalance
    ! integer(ccs_int) :: seed
    ! integer(ccs_int) :: mode
    ! integer(ccs_int) :: suppress
    ! integer(ccs_int) :: edgecuts
    integer(ccs_int) :: local_part_size
    integer(ccs_int) :: irank
    ! integer(ccs_int) :: isize
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

    ! Values hardcoded for now
    imbalance = 0.03  ! Desired balance - 0.03 = 3%
    seed = 2022       ! "Random" seed
    mode = 4          ! FASTSOCIAL
    suppress = 0      ! Do not suppress the output
    edgecuts = -1     ! XXX: silence unused variable warning

    allocate (tmp_partition(mesh%topo%global_num_cells)) ! Temporary partition array
    tmp_partition = 0

    irank = par_env%proc_id ! Current rank
    ! isize = par_env%num_procs
    num_procs = par_env%num_procs

    ! ParHIP needs 0-indexing - shift array contents by -1
    ! mesh%topo%vtxdist = mesh%topo%vtxdist - 1
    ! mesh%topo%xadj = mesh%topo%xadj - 1
    ! mesh%topo%adjncy = mesh%topo%adjncy - 1

    vtxdist = mesh%topo%vtxdist - 1
    xadj = mesh%topo%xadj - 1
    adjncy = mesh%topo%adjncy - 1
    
    ! Set weights to 1
    ! mesh%topo%adjwgt = 1
    ! mesh%topo%vwgt = 1
    adjwgt = 1
    vwgt = 1

    ! Number of elements in local partition array
    ! Needed for gathering loca partitions into global partition array
    local_part_size = size(mesh%topo%local_partition)

    allocate(local_partition(local_part_size))

    ! Partitioning an unweighted graph
    select type (par_env)
    type is (parallel_environment_mpi)

      ! call partition_parhipkway(mesh%topo%vtxdist, mesh%topo%xadj, mesh%topo%adjncy, &
      !                           mesh%topo%vwgt, mesh%topo%adjwgt, &
      !                           par_env%num_procs, imbalance, suppress, &
      !                           seed, mode, edgecuts, mesh%topo%local_partition, par_env%comm)

      call partition_parhipkway(vtxdist, xadj, adjncy, vwgt, adjwgt, &
                                num_procs, imbalance, suppress, &
                                seed, mode, edgecuts, local_partition, comm)

      mesh%topo%local_partition(:) = local_partition(:)

      do i = 1, local_part_size
        tmp_partition(i + mesh%topo%vtxdist(irank + 1)) = mesh%topo%local_partition(i)
      end do

      call MPI_AllReduce(tmp_partition, mesh%topo%global_partition, mesh%topo%global_num_cells, &
                         MPI_LONG, MPI_SUM, par_env%comm, ierr)

    class default
      print *, "ERROR: Unknown parallel environment!"
    end select

    call dprint("Number of edgecuts: " // str(edgecuts))

    ! Return to 1-indexing by adding 1
    ! mesh%topo%vtxdist = mesh%topo%vtxdist + 1
    ! mesh%topo%xadj = mesh%topo%xadj + 1
    ! mesh%topo%adjncy = mesh%topo%adjncy + 1

    deallocate (tmp_partition)

  end subroutine

end module parhip_interfaces
