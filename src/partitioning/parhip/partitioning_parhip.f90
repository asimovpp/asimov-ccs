submodule (partitioning) partitioning_parhip
#include "ccs_macros.inc"

  use kinds, only: ccs_int, ccs_real
  use utils, only : str, debug_print
  use parallel_types_mpi, only: parallel_environment_mpi

  implicit none

contains

  !v Partition the mesh
  !
  ! Use ParHIP library to compute a k-way vertex separator given a k-way partition of the graph.
  ! The graph can be weighted or unweighted.
  module subroutine partition_kway(par_env, topo)

    use mpi

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(topology), target, intent(inout) :: topo                           !< The topology for which to compute the parition

    ! Local variables
    integer(ccs_long), dimension(:), allocatable :: tmp_partition
    real(ccs_real) :: imbalance
    integer(ccs_int) :: seed
    integer(ccs_int) :: mode
    integer(ccs_int) :: suppress
    integer(ccs_int) :: edgecuts
    integer(ccs_int) :: local_part_size
    integer(ccs_int) :: irank
    integer(ccs_int) :: ierr
    integer(ccs_int) :: i

    ! Values hardcoded for now
    imbalance = 0.03  ! Desired balance - 0.03 = 3% 
    seed = 2022       ! "Random" seed
    mode = 4          ! FASTSOCIAL
    suppress = 0      ! Do not suppress the output
    edgecuts = -1     ! XXX: silence unused variable warning

    allocate(tmp_partition(topo%global_num_cells)) ! Temporary partition array
    tmp_partition = 0

    irank = par_env%proc_id ! Current rank

    ! ParHIP needs 0-indexing - shift array contents by -1
    topo%vtxdist = topo%vtxdist - 1
    topo%xadj = topo%xadj - 1
    topo%adjncy = topo%adjncy - 1

    ! Set weights to 1
    topo%adjwgt = 1
    topo%vwgt = 1

    ! Number of elements in local partition array
    ! Needed for gathering loca partitions into global partition array
    local_part_size = size(topo%local_partition)

    ! Partitioning an unweighted graph
    select type(par_env)
    type is (parallel_environment_mpi)

    call partition_parhipkway(topo%vtxdist, topo%xadj, topo%adjncy, &
                              topo%vwgt, topo%adjwgt, & 
                              par_env%num_procs, imbalance, suppress, &
                              seed, mode, edgecuts, topo%local_partition, par_env%comm)

    do i=1,local_part_size
      tmp_partition(i+topo%vtxdist(irank+1)) = topo%local_partition(i)
    end do

    call MPI_AllReduce(tmp_partition,topo%global_partition,topo%global_num_cells,MPI_LONG,MPI_SUM,par_env%comm,ierr)

    class default
      print*, "ERROR: Unknown parallel environment!"
    end select    

    call dprint ("Number of edgecuts: "//str(edgecuts))

    ! Return to 1-indexing by adding 1
    topo%vtxdist = topo%vtxdist + 1
    topo%xadj = topo%xadj + 1
    topo%adjncy = topo%adjncy + 1

    deallocate(tmp_partition)
    
  end subroutine partition_kway

  !v Compute the input arrays for the partitioner
  !
  ! Using the topology object, compute the input arrays for the ParHIP partitioner
  ! Input arrays for the partitioner are: vtxdist, xadj and adjncy
  module subroutine compute_partitioner_input(par_env, topo)

    use iso_fortran_env, only: int32

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(topology), target, intent(inout) :: topo                           !< The topology for which to compute the parition

    ! Local variables
    integer(ccs_int), dimension(:,:), allocatable :: tmp_int2d ! Temporary 2D integer array

    integer(ccs_int) :: i, j, k
    integer(ccs_int) :: irank ! MPI rank ID
    integer(ccs_int) :: isize ! Size of MPI world
    integer(ccs_int) :: start_index 
    integer(ccs_int) :: end_index  
    integer(ccs_int) :: face_nb1  
    integer(ccs_int) :: face_nb2
    integer(ccs_int) :: local_index
    integer(ccs_int) :: num_connections
 
    irank = par_env%proc_id
    isize = par_env%num_procs
  
    ! Create and populate the vtxdist array based on the total number of cells
    ! and the total number of ranks in the parallel environment
    allocate(topo%vtxdist(isize + 1)) ! vtxdist array is of size num_procs + 1 on all ranks

    topo%vtxdist(1) = 1                                 ! First element is 1
    topo%vtxdist(isize + 1) = topo%global_num_cells + 1 ! Last element is total number of cells + 1

    ! Divide the total number of cells by the world size to
    ! compute the chunk sizes
    k = int(real(topo%global_num_cells) / isize)
    j = 1

    do i = 1, isize
      topo%vtxdist(i) = j
      j = j + k
    enddo

    start_index = int(topo%vtxdist(irank + 1), int32)
    end_index = int(topo%vtxdist(irank + 2), int32) - 1
  
    ! Allocate adjacency array xadj based on vtxdist
    allocate(topo%xadj(topo%vtxdist(irank + 2) - topo%vtxdist(irank + 1) + 1)) 

    ! Allocated temporary 2D integer work array and initialise to 0
    allocate(tmp_int2d(topo%vtxdist(irank + 2) - topo%vtxdist(irank + 1), topo%max_faces + 1))
    tmp_int2d = 0

    ! Allocate global boundaries array
    allocate(topo%global_boundaries(topo%global_num_cells))

    ! All ranks loop over all the faces
    do i=1,topo%global_num_faces

      face_nb1 = topo%face_edge_end1(i)
      face_nb2 = topo%face_edge_end2(i)

      ! If face neighbour 1 is local to the current rank
      ! and face neighbour 2 is not 0
      if(face_nb1 .ge. start_index .and. face_nb1 .le. end_index .and. face_nb2 .ne. 0) then
         local_index = face_nb1 - start_index + 1                 ! Local cell index
         k = tmp_int2d(local_index, topo%max_faces + 1) + 1  ! Increment number of faces for this cell
         tmp_int2d(local_index, k) = face_nb2                       ! Store global index of neighbour cell
         tmp_int2d(local_index, topo%max_faces + 1) = k      ! Store number of faces for this cell
      endif

      ! If face neighbour 2 is local to the current rank
      ! and face neighbour 1 is not 0
      if(face_nb2 .ge. start_index .and. face_nb2 .le. end_index .and. face_nb1  .ne. 0) then
         local_index = face_nb2 - start_index + 1                 ! Local cell index
         k = tmp_int2d(local_index, topo%max_faces + 1) + 1  ! Increment number of faces for this cell
         tmp_int2d(local_index, k) = face_nb1                       ! Store global index of neighbour cell
         tmp_int2d(local_index, topo%max_faces + 1) = k      ! Store number of faces for this cell
      endif

      ! If face neighbour 2 is 0 we have a boundary face
      if(face_nb2 .eq. 0) then
        topo%global_boundaries(face_nb1) = topo%global_boundaries(face_nb1) + 1
      end if

    enddo 

    num_connections = sum(tmp_int2d(:, topo%max_faces+1))

    ! Allocate adjncy array based on the number of computed connections
    allocate(topo%adjncy(num_connections))
    ! Allocate local partition array
    allocate(topo%local_partition(topo%vtxdist(irank+2)-topo%vtxdist(irank+1)))
    ! Allocate global partition array
    allocate(topo%global_partition(topo%global_num_cells))

    local_index = 1

    do i=1,end_index-start_index+1  ! Loop over local cells
      
      topo%xadj(i) = local_index                          ! Pointer to start of current
       
      do j=1,tmp_int2d(i, topo%max_faces+1)               ! Loop over number of faces
        topo%adjncy(local_index + j - 1) = tmp_int2d(i,j) ! Store global IDs of neighbour cells
      end do

       local_index = local_index + tmp_int2d(i,topo%max_faces+1)
       topo%xadj(i+1) = local_index

    end do

    ! Allocate weight arrays
    allocate(topo%adjwgt(num_connections))
    ! Allocate local partition array
    allocate(topo%vwgt(topo%vtxdist(irank+2)-topo%vtxdist(irank+1)))    

    deallocate(tmp_int2d)

  end subroutine compute_partitioner_input

end submodule