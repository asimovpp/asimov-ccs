submodule(partitioning) partitioning_parhip
#include "ccs_macros.inc"

  use kinds, only: ccs_int, ccs_real, ccs_long
  use utils, only: str, debug_print
  use parallel_types_mpi, only: parallel_environment_mpi
  use meshing, only: set_local_num_cells, get_local_num_cells

  implicit none

  interface 
    subroutine partition_parhipkway(vtxdist, xadj, adjncy, vwgt, adjwgt, &
                                  num_procs, imbalance, suppress, &
                                  seed, mode, edgecuts, local_partition, comm) bind(c)
      use iso_c_binding

      integer(c_long) :: vtxdist(*)
      integer(c_long) :: xadj(*)
      integer(c_long) :: adjncy(*)
      integer(c_long) :: vwgt(*)
      integer(c_long) :: adjwgt(*)
      integer(c_int) :: num_procs
      real(c_double) :: imbalance
      integer(c_int) :: suppress
      integer(c_int) :: seed
      integer(c_int) :: mode
      integer(c_int) :: edgecuts
      integer(c_long) :: local_partition(*)
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
    integer(ccs_long), dimension(:), allocatable :: tmp_partition
    integer(ccs_int) :: local_part_size
    integer(ccs_int) :: irank
    integer(ccs_int) :: ierr
    integer(ccs_int) :: i

    integer(c_long), allocatable :: vtxdist(:)
    integer(c_long), allocatable :: xadj(:)
    integer(c_long), allocatable :: adjncy(:)
    integer(c_long), allocatable :: vwgt(:)
    integer(c_long), allocatable :: adjwgt(:)
    integer(c_int) :: num_procs
    real(c_double) :: imbalance
    integer(c_int) :: suppress
    integer(c_int) :: seed
    integer(c_int) :: mode
    integer(c_int) :: edgecuts
    integer(c_long), allocatable :: local_partition(:)
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

    allocate(local_partition(local_part_size))

    ! Partitioning an unweighted graph
    select type (par_env)
    type is (parallel_environment_mpi)

      comm = par_env%comm

      call partition_parhipkway(vtxdist, xadj, adjncy, vwgt, adjwgt, &
                                num_procs, imbalance, suppress, &
                                seed, mode, edgecuts, local_partition, comm)

      mesh%topo%local_partition(:) = local_partition(:)

      do i = 1, local_part_size
        tmp_partition(i + vtxdist(irank+1)) = mesh%topo%local_partition(i)
      end do

      call MPI_AllReduce(tmp_partition, mesh%topo%global_partition, mesh%topo%global_num_cells, &
                         MPI_LONG, MPI_SUM, par_env%comm, ierr)

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
