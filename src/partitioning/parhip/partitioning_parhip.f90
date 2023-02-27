submodule(partitioning) partitioning_parhip
#include "ccs_macros.inc"

  use kinds, only: ccs_int, ccs_real
  use utils, only: str, debug_print
  use parallel_types_mpi, only: parallel_environment_mpi
  use meshing, only: set_local_num_cells, get_local_num_cells

  implicit none

contains

  !v Partition the mesh
  !
  ! Use ParHIP library to compute a k-way vertex separator given a k-way partition of the graph.
  ! The graph can be weighted or unweighted.
  module subroutine partition_kway(par_env, mesh)

    use mpi

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the parition

    ! Local variables
    integer(ccs_long), dimension(:), allocatable :: tmp_partition
    real(ccs_real) :: imbalance
    integer(ccs_int) :: seed
    integer(ccs_int) :: mode
    integer(ccs_int) :: suppress
    integer(ccs_int) :: edgecuts
    integer(ccs_int) :: local_part_size
    integer(ccs_int) :: irank, isize
    integer(ccs_int) :: ierr
    integer(ccs_int) :: i
    
    ! Values hardcoded for now
    imbalance = 0.03  ! Desired balance - 0.03 = 3%
    seed = 2022       ! "Random" seed
    mode = 4          ! FASTSOCIAL
    suppress = 0      ! Do not suppress the output
    edgecuts = -1     ! XXX: silence unused variable warning

    allocate (tmp_partition(mesh%topo%global_num_cells)) ! Temporary partition array
    tmp_partition = 0

    irank = par_env%proc_id ! Current rank
    isize = par_env%num_procs

    ! ParHIP needs 0-indexing - shift array contents by -1
    mesh%topo%vtxdist = mesh%topo%vtxdist - 1
    mesh%topo%xadj = mesh%topo%xadj - 1
    mesh%topo%adjncy = mesh%topo%adjncy - 1

    ! Set weights to 1
    mesh%topo%adjwgt = 1
    mesh%topo%vwgt = 1

    ! Number of elements in local partition array
    ! Needed for gathering loca partitions into global partition array
    local_part_size = size(mesh%topo%local_partition)

    ! Partitioning an unweighted graph
    select type (par_env)
    type is (parallel_environment_mpi)

      call partition_parhipkway(mesh%topo%vtxdist, mesh%topo%xadj, mesh%topo%adjncy, &
                                mesh%topo%vwgt, mesh%topo%adjwgt, &
                                par_env%num_procs, imbalance, suppress, &
                                seed, mode, edgecuts, mesh%topo%local_partition, par_env%comm)

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
    mesh%topo%vtxdist = mesh%topo%vtxdist + 1
    mesh%topo%xadj = mesh%topo%xadj + 1
    mesh%topo%adjncy = mesh%topo%adjncy + 1

    deallocate (tmp_partition)

  end subroutine partition_kway

  !v Compute the input arrays for the partitioner
  !
  ! Using the topology object, compute the input arrays for the ParHIP partitioner
  ! Input arrays for the partitioner are: vtxdist, xadj and adjncy
  module subroutine compute_partitioner_input(par_env, mesh)

    use iso_fortran_env, only: int32

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the parition

    ! Local variables
    integer(ccs_int), dimension(:, :), allocatable :: tmp_int2d ! Temporary 2D integer array

    integer(ccs_int) :: i, j, k
    integer(ccs_int) :: irank ! MPI rank ID
    integer(ccs_int) :: isize ! Size of MPI world
    integer(ccs_int) :: start_index
    integer(ccs_int) :: end_index
    integer(ccs_int) :: face_nb1
    integer(ccs_int) :: face_nb2
    integer(ccs_int) :: local_index
    integer(ccs_int) :: num_connections
    integer(ccs_int) :: local_num_cells

    irank = par_env%proc_id
    isize = par_env%num_procs

    start_index = int(mesh%topo%vtxdist(irank + 1), int32)
    end_index = int(mesh%topo%vtxdist(irank + 2), int32) - 1

    ! Allocate global partition array
    allocate (mesh%topo%global_partition(mesh%topo%global_num_cells))

    ! Initial global partition
    do i = 1, size(mesh%topo%vtxdist) - 1
      j = i - 1
      mesh%topo%global_partition(mesh%topo%vtxdist(i):mesh%topo%vtxdist(i + 1) - 1) = j
    end do

    ! Count the number of local cells per rank
    local_num_cells = count(mesh%topo%global_partition == irank)
    call set_local_num_cells(local_num_cells, mesh)
    call get_local_num_cells(mesh, local_num_cells) ! Ensure using true value
    call dprint("Initial number of local cells: " // str(local_num_cells))

    ! Allocate adjacency index array xadj based on vtxdist
    allocate (mesh%topo%xadj(mesh%topo%vtxdist(irank + 2) - mesh%topo%vtxdist(irank + 1) + 1))

    ! Allocate temporary 2D integer work array and initialise to 0
    allocate (tmp_int2d(mesh%topo%vtxdist(irank + 2) - mesh%topo%vtxdist(irank + 1), mesh%topo%max_faces + 1))
    tmp_int2d = 0

    ! Allocate global boundaries array
    allocate (mesh%topo%global_boundaries(mesh%topo%global_num_cells))

    ! All ranks loop over all the faces
    do i = 1, mesh%topo%global_num_faces

      face_nb1 = mesh%topo%face_cell1(i)
      face_nb2 = mesh%topo%face_cell2(i)

      ! If face neighbour 1 is local to the current rank
      ! and face neighbour 2 is not 0
      if (face_nb1 .ge. start_index .and. face_nb1 .le. end_index .and. face_nb2 .ne. 0) then
        local_index = face_nb1 - start_index + 1                 ! Local cell index
        k = tmp_int2d(local_index, mesh%topo%max_faces + 1) + 1  ! Increment number of faces for this cell
        tmp_int2d(local_index, k) = face_nb2                       ! Store global index of neighbour cell
        tmp_int2d(local_index, mesh%topo%max_faces + 1) = k      ! Store number of faces for this cell
      end if

      ! If face neighbour 2 is local to the current rank
      ! and face neighbour 1 is not 0
      if (face_nb2 .ge. start_index .and. face_nb2 .le. end_index .and. face_nb1 .ne. 0) then
        local_index = face_nb2 - start_index + 1                 ! Local cell index
        k = tmp_int2d(local_index, mesh%topo%max_faces + 1) + 1  ! Increment number of faces for this cell
        tmp_int2d(local_index, k) = face_nb1                       ! Store global index of neighbour cell
        tmp_int2d(local_index, mesh%topo%max_faces + 1) = k      ! Store number of faces for this cell
      end if

      ! If face neighbour 2 is 0 we have a boundary face
      if (face_nb2 .eq. 0) then
        mesh%topo%global_boundaries(face_nb1) = mesh%topo%global_boundaries(face_nb1) + 1
      end if

    end do

    num_connections = sum(tmp_int2d(:, mesh%topo%max_faces + 1))
    call dprint("Initial number of connections: " // str(num_connections))

    ! Allocate adjncy array based on the number of computed connections
    allocate (mesh%topo%adjncy(num_connections))
    ! Allocate local partition array
    allocate (mesh%topo%local_partition(mesh%topo%vtxdist(irank + 2) - mesh%topo%vtxdist(irank + 1)))

    local_index = 1

    do i = 1, end_index - start_index + 1  ! Loop over local cells

      mesh%topo%xadj(i) = local_index                          ! Pointer to start of current

      do j = 1, tmp_int2d(i, mesh%topo%max_faces + 1)               ! Loop over number of faces
        mesh%topo%adjncy(local_index + j - 1) = tmp_int2d(i, j) ! Store global IDs of neighbour cells
        if (mesh%topo%global_partition(tmp_int2d(i, j)) /= irank) then
          mesh%topo%halo_num_cells = mesh%topo%halo_num_cells + 1
        end if
      end do

      local_index = local_index + tmp_int2d(i, mesh%topo%max_faces + 1)
      mesh%topo%xadj(i + 1) = local_index

    end do

    call dprint("Initial number of halo cells: " // str(mesh%topo%halo_num_cells))

    mesh%topo%total_num_cells = local_num_cells + mesh%topo%halo_num_cells

    call dprint("Total number of cells (local + halo): " // str(mesh%topo%total_num_cells))

    ! Allocate weight arrays
    allocate (mesh%topo%adjwgt(num_connections))
    ! Allocate local partition array
    allocate (mesh%topo%vwgt(mesh%topo%vtxdist(irank + 2) - mesh%topo%vtxdist(irank + 1)))

    deallocate (tmp_int2d)

  end subroutine compute_partitioner_input

end submodule
