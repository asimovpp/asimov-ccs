submodule (partitioning) partitioning_common
#include "ccs_macros.inc"

  use kinds, only: ccs_int
  use utils, only : str, debug_print
  use parallel_types_mpi, only: parallel_environment_mpi

implicit none

  contains

  !v Read the topology data from an input (HDF5) file
  ! This subroutine assumes the following names are used in the file:
  ! "ncel" - the total number of cells
  ! "nfac" - the total number of faces
  ! "maxfaces" - the maximum number of faces per cell
  ! "/face/cell1" and "/face/cell2" - the arrays the face edge data
  module subroutine read_topology(par_env, case_name, topo)

    use constants, only: geoext, adiosconfig
    use io, only: initialise_io, cleanup_io, configure_io, &
                  open_file, close_file, &
                  read_scalar, read_array
    use types, only: io_environment, io_process
    use parallel, only: read_command_line_arguments

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    character(len=:), allocatable :: case_name                              !< The name of the case that is computed
    type(topology), target, intent(inout) :: topo                           !< The topology for which to compute the parition

    ! Local variables
    character(len=:), allocatable :: geo_file    ! Geo file name
    character(len=:), allocatable :: adios2_file ! ADIOS2 config file name

    class(io_environment), allocatable :: io_env
    class(io_process), allocatable :: geo_reader

    integer(ccs_long), dimension(1) :: sel_start
    integer(ccs_long), dimension(1) :: sel_count

    integer(ccs_long), dimension(2) :: sel2_start
    integer(ccs_long), dimension(2) :: sel2_count

    geo_file = case_name//geoext
    adios2_file = case_name//adiosconfig

    call initialise_io(par_env, adios2_file, io_env)
    call configure_io(io_env, "geo_reader", geo_reader)  
  
    call open_file(geo_file, "read", geo_reader)
  
    ! Read attribute "ncel" - the total number of cells
    call read_scalar(geo_reader, "ncel", topo%global_num_cells)
    ! Read attribute "nfac" - the total number of faces
    call read_scalar(geo_reader, "nfac", topo%global_num_faces)
    ! Read attribute "maxfaces" - the maximum number of faces per cell
    call read_scalar(geo_reader, "maxfaces", topo%max_faces)

    allocate(topo%face_cell1(topo%global_num_faces))
    allocate(topo%face_cell2(topo%global_num_faces))
    allocate(topo%global_face_indices(topo%max_faces, topo%global_num_cells))

    sel_start(1) = 0 ! Global index to start reading from
    sel_count(1) = topo%global_num_faces ! How many elements to read in total

    ! Read arrays face/cell1 and face/cell2
    call read_array(geo_reader, "/face/cell1", sel_start, sel_count, topo%face_cell1)
    call read_array(geo_reader, "/face/cell2", sel_start, sel_count, topo%face_cell2)

    sel2_start = 0
    sel2_count(1) = topo%max_faces! topo%global_num_cells
    sel2_count(2) = topo%global_num_cells

    call read_array(geo_reader, "/cell/cface", sel2_start, sel2_count, topo%global_face_indices)

    ! Close the file and ADIOS2 engine
    call close_file(geo_reader)

    ! Finalise the ADIOS2 IO environment
    call cleanup_io(io_env)

  end subroutine read_topology

  ! Compute the new topology connectivity after partitioning
  module subroutine compute_connectivity(par_env, topo)

    use mpi
    use iso_fortran_env, only: int32

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(topology), target, intent(inout) :: topo                           !< The topology for which to compute the parition

    ! Local variables
    integer(ccs_int), dimension(:,:), allocatable :: tmp_int2d ! Temporary 2D integer array
    integer(ccs_int) :: irank ! MPI rank ID
    integer(ccs_int) :: isize ! Size of MPI world
    integer(ccs_int) :: i, j, k, n
    integer(ccs_int) :: start_index 
    integer(ccs_int) :: end_index  
    integer(ccs_int) :: face_nb1  
    integer(ccs_int) :: face_nb2
    integer(ccs_int) :: local_index
    integer(ccs_int) :: num_connections
    integer(ccs_int) :: current, previous
 
    irank = par_env%proc_id
    isize = par_env%num_procs

    ! Count the new number of local cells per rank
    topo%local_num_cells = count(topo%global_partition == irank)
    call dprint ("Number of local cells after partitioning: "//str(topo%local_num_cells))

    topo%vtxdist(1) = 1

    ! Recompute vtxdist array based on the new partition
    do i = 2, isize + 1
      topo%vtxdist(i) = count(topo%global_partition == i - 2) + topo%vtxdist(i - 1)
    end do

    if(irank == 0) then
      do i = 1, isize + 1
        call dprint("new vtxdist("//str(i)//"): "//str(int(topo%vtxdist(i))))
      end do    
    end if

    ! Deallocate old xadj array
    if (allocated(topo%xadj)) then
      deallocate(topo%xadj)
    end if

    ! Allocate new adjacency index array xadj based on new vtxdist
    allocate(topo%xadj(topo%vtxdist(irank + 2) - topo%vtxdist(irank + 1) + 1)) 
  
    start_index = int(topo%vtxdist(irank + 1), int32)
    end_index = int(topo%vtxdist(irank + 2), int32) - 1

    ! Allocate array to hold number of neighbours for local cells
    allocate(topo%num_nb(topo%local_num_cells))
    topo%num_nb(:) = topo%max_faces ! NOTE: This assumes a mesh of uniform cell types

    call dprint("Number of neighbours "//str(int(topo%num_nb(1))))

    if(allocated(topo%global_boundaries) .eqv. .false.) then
      allocate(topo%global_boundaries(topo%global_num_cells))
    end if

    ! Reset global_boundaries array
    topo%global_boundaries = 0

   ! Allocate temporary 2D integer work array and initialise to 0
    allocate(tmp_int2d(topo%vtxdist(irank + 2) - topo%vtxdist(irank + 1), topo%max_faces + 1))
    tmp_int2d = 0

    ! All ranks loop over all the faces again
    do i=1,topo%global_num_faces

      face_nb1 = topo%face_cell1(i)
      face_nb2 = topo%face_cell2(i)

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

    ! New number of local connections
    num_connections = sum(tmp_int2d(:, topo%max_faces+1))
    call dprint ("Number of connections after partitioning: "//str(num_connections))

    ! Allocate new adjncy array based on the new number of computed connections
    if (allocated(topo%adjncy)) then
      deallocate(topo%adjncy)
    end if
    
    allocate(topo%adjncy(num_connections))

    local_index = 1

    do i=1,end_index-start_index+1  ! Loop over local cells
      
      topo%xadj(i) = local_index                          ! Pointer to start of current
       
      do j=1,tmp_int2d(i, topo%max_faces+1)               ! Loop over number of faces
        topo%adjncy(local_index + j - 1) = tmp_int2d(i,j) ! Store global IDs of neighbour cells
        if(topo%global_partition(tmp_int2d(i,j)) /= irank) then
          topo%halo_num_cells = topo%halo_num_cells + 1   ! Count the number of halo cells
        end if
      end do

       local_index = local_index + tmp_int2d(i,topo%max_faces+1)
       topo%xadj(i+1) = local_index

    end do

    call dprint ("Number of halo cells after partitioning: "//str(topo%halo_num_cells))

    topo%total_num_cells = topo%local_num_cells + topo%halo_num_cells

    call dprint ("Total number of cells (local + halo) after partitioning: "//str(topo%total_num_cells))

    ! Allocate and then compute global indices
    if(allocated(topo%global_indices) .eqv. .false.) then
      allocate(topo%global_indices(topo%total_num_cells))
    end if

    ! First the indices of the local cells
    k = 1
    do i = 1, topo%global_num_cells
      if(topo%global_partition(i) == irank) then
        topo%global_indices(k) = i
        call dprint("Index of local cell: topo%global_indices("//str(k)//") = "//str(i))
        k = k + 1
      end if 
    end do
    ! Then the indices of the halo cells
    k = topo%local_num_cells + 1
    do i = 1, size(topo%adjncy)
      if (topo%global_partition(topo%adjncy(i)) /= irank) then
        call dprint("Index of halo cell: topo%global_indices("//str(k)//") = "//str(i))
        topo%global_indices(k) = topo%adjncy(i)
        k = k + 1
      end if
    end do

    ! Now compute the face indices for the local vector
    allocate(topo%face_indices(topo%max_faces, topo%local_num_cells))
    do i = 1, topo%local_num_cells
      k = topo%global_indices(i)
      do j = 1, topo%max_faces
        topo%face_indices(j,i) = topo%global_face_indices(j,k) ! XXX: currently using global faces numbers
      end do
    end do

    do i = 1, size(topo%global_indices) 
      call dprint("Global index "//str(i)//": "//str(int(topo%global_indices(i))))
    end do

    deallocate(topo%global_face_indices) ! No longer needed now we have the local face indices

    topo%num_faces = 0
    previous = 0

    ! Next, compute the number of local faces, which  
    ! equals the maximum value in the face_indices array
    topo%num_faces = maxval(topo%face_indices)

    call dprint("Number of local faces "//str(topo%num_faces))

    if(irank == 0) print*,"Global boundaries: ", topo%global_boundaries(1:10)

    ! Finally, allocate and compute the neighbour indices
    if (allocated(topo%nb_indices)) then
      deallocate(topo%nb_indices)
    end if
    allocate(topo%nb_indices(topo%max_faces, topo%local_num_cells))

    do i = 1, size(topo%xadj) - 1
      n = topo%xadj(i+1) - topo%xadj(i)
      do j = 1, n
        topo%nb_indices(j, i) = topo%adjncy(j+topo%xadj(i)-1)
      end do 
      do k = n + 1, topo%max_faces
        topo%nb_indices(k, i) = 0 ! Set boundaries to 0 - will be updated with correct BCs
      end do
    end do

    if (irank == 0) print*, "Neighbour indices ", topo%nb_indices

  end subroutine compute_connectivity

end submodule