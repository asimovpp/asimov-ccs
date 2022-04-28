submodule (partitioning) partitioning_parhip

  use kinds, only: ccs_int, ccs_real

  implicit none

contains

  !v Partition the mesh
  !
  ! Use ParHIP library to compute a k-way vertex separator given a k-way partition of the graph.
  ! The graph can be weighted or unweighted.
  module subroutine partition_kway(par_env, topo, partition_vector)

    use iso_c_binding

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(topology), target, intent(inout) :: topo                           !< The topology for which to compute the parition
    integer(ccs_long), allocatable, intent(out) :: partition_vector         !< The vector with the partition 

    ! Local variables
    real(ccs_real) :: imbalance
    integer(ccs_int) :: seed
    integer(ccs_int) :: mode
    integer(ccs_int) :: suppress
    integer(ccs_int) :: edgecut
    type(c_ptr) :: vwgt     ! NULL pointer to be used if no weights are attached to vertices
    type(c_ptr) :: adjwgt   ! NULL pointer to be used if no weights are attached to edges

    ! Values hardcoded for now
    imbalance = 0.03  ! Desired balance - 0.03 = 3% 
    seed = 2022       ! "Random" seed
    mode = 3          ! FASTSOCIAL
    suppress = 0      ! Do not suppress the output
    edgecut = -1      ! XXX: silence unused variable warning

    ! ParHIP needs 0-indexing
    topo%vtxdist = topo%vtxdist -1
    topo%xadj = topo%xadj - 1
    topo%adjncy = topo%adjncy - 1

    ! If topo%vwgt and topo%adjwgt are not allocated, no weights have been set
    ! and we are dealing with an unweighted graph
    if(allocated(topo%vwgt) .eqv. .false. .and. allocated(topo%adjwgt) .eqv. .false.) then

      ! NULL pointers
      vwgt = c_null_ptr
      adjwgt = c_null_ptr
   
      ! Partitioning an unweighted graph
      ! call partition_parhipkway(topo%vtxdist, topo%xadj, topo%adjncy, vwgt, adjwgt, & 
      !                             par_env%num_procs, imbalance, suppress, seed, &
      !                             mode, edgecut, partition_vector, par_env%comm)
      
    ! If topo%vwgt and topo%adjwgt are allocated we are dealing with an weighted graph
    else 

      ! Partitioning a graph with weights on vertices and edges
      ! call partition_parhipkway(topo%vtxdist, topo%xadj, topo%adjncy, topo%vwgt, topo%adjwgt, & 
      !                             par_env%num_procs, imbalance, suppress, seed, &
      !                             mode, edgecut, partition_vector, par_env%comm)

    end if

    ! Return to 1-indexing. For safety only right now, later on these arrays
    ! will be deallocated as they are not used beyond this point
    topo%vtxdist = topo%vtxdist -1
    topo%xadj = topo%xadj - 1
    topo%adjncy = topo%adjncy - 1

    ! Make sure partition vector starts from 1
    partition_vector = partition_vector + 1

  end subroutine partition_kway

  !v Compute the input arrays for the partitioner
  !
  ! Using the topology object, compute the input arrays for the ParHIP partitioner
  module subroutine compute_partitioner_input(par_env, topo)

    use constants, only: geoext, adiosconfig
    use io, only: initialise_io, cleanup_io, configure_io, &
                  open_file, close_file, &
                  read_scalar, read_array
    use types, only: io_environment, io_process
    use parallel, only: read_command_line_arguments

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(topology), target, intent(inout) :: topo                           !< The topology for which to compute the parition

    ! Local variables
    character(len=:), allocatable :: case_name   ! Case name
    character(len=:), allocatable :: geo_file    ! Geo file name
    character(len=:), allocatable :: adios2_file ! ADIOS2 config file name

    class(io_environment), allocatable :: io_env
    class(io_process), allocatable :: geo_reader

    integer(ccs_int) :: i, j, k
    integer(ccs_int) :: irank ! MPI rank ID
    integer(ccs_int) :: isize ! Size of MPI world
    integer(ccs_int) :: start_index 
    integer(ccs_int) :: end_index  

    integer(ccs_int) :: num_cells ! Total number of cells
    integer(ccs_int) :: num_faces ! Total number of faces
  
    irank = par_env%proc_id
    isize = par_env%num_procs
  
    call read_command_line_arguments(par_env, case_name=case_name)
  
    geo_file = case_name//geoext
    adios2_file = case_name//adiosconfig     

    call initialise_io(par_env, adios2_file, io_env)
    call configure_io(io_env, "geo_reader", geo_reader)  
  
    call open_file(geo_file, "read", geo_reader)
  
    ! Read attribute "ncel" - the total number of cells
    call read_scalar(geo_reader, "ncel", num_cells)
    ! Read attribute "nfac" - the total number of faces
    call read_scalar(geo_reader, "nfac", num_faces)

    ! Create and populate the vtxdist array based on the total number of cells
    ! and the total number of ranks in the parallel environment
    allocate(topo%vtxdist(isize + 1)) ! vtxdist array is of size num_procs + 1 on all ranks

    topo%vtxdist(1) = 1                     ! First element is 1
    topo%vtxdist(isize + 1) = num_cells + 1 ! Last element is total number of cells + 1

    ! Divide the total number of cells by the world size to
    ! compute the chunk sizes
    k = int(real(num_cells) / isize)
    j = 1

    do i = 1, isize
      topo%vtxdist(i) = j
      j = j + k
    enddo

    start_index = topo%vtxdist(irank+1)
    end_index = topo%vtxdist(irank+2)-1
  
    ! Create xadj array based on the content computed for vtxdist
    allocate(topo%xadj(topo%vtxdist(irank+2)-topo%vtxdist(irank+1)+1)) 
  
    ! TODO - REMOVE DEALLOCATES - ARRAYS ARE NEEDED OUTSIDE THIS SUBROUTINE
    deallocate(topo%vtxdist)
    deallocate(topo%xadj)

  end subroutine compute_partitioner_input

end submodule