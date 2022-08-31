submodule (partitioning) partitioning_common
#include "ccs_macros.inc"

  use kinds, only: ccs_int
  use utils, only : str, debug_print
  use parallel_types_mpi, only: parallel_environment_mpi

implicit none

  contains


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
    integer(ccs_int) :: i
    integer(ccs_int) :: start_index 
    integer(ccs_int) :: end_index  
    integer(ccs_int) :: face_nb1  
    integer(ccs_int) :: face_nb2
    integer(ccs_int) :: num_connections
 
    irank = par_env%proc_id
    isize = par_env%num_procs

    ! Count the new number of local cells per rank
    topo%local_num_cells = count(topo%global_partition == irank)
    call dprint ("Number of local cells after partitioning: "//str(topo%local_num_cells))

    ! Get global indices of local cells
    call compute_connectivity_get_local_cells(par_env, topo)

    ! Recompute vtxdist array based on the new partition
    topo%vtxdist(1) = 1
    do i = 2, isize + 1
      topo%vtxdist(i) = count(topo%global_partition == (i - 2)) + topo%vtxdist(i - 1)
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

    if(allocated(topo%global_boundaries) .eqv. .false.) then
      allocate(topo%global_boundaries(topo%global_num_cells))
    end if

    ! Reset global_boundaries array
    topo%global_boundaries = 0

    ! Allocate temporary 2D integer work array and initialise to 0
    allocate(tmp_int2d(topo%vtxdist(irank + 2) - topo%vtxdist(irank + 1), topo%max_faces + 1))
    tmp_int2d = 0

    start_index = int(topo%vtxdist(irank + 1), int32)
    end_index = int(topo%vtxdist(irank + 2), int32) - 1

    ! Allocate array to hold number of neighbours for local cells
    if (allocated(topo%num_nb)) then
      deallocate(topo%num_nb)
    end if
    allocate(topo%num_nb(topo%local_num_cells))
    
    ! All ranks loop over all the faces again
    do i=1,topo%global_num_faces

      face_nb1 = topo%face_cell1(i)
      face_nb2 = topo%face_cell2(i)

      ! If face neighbour 1 is local to the current rank
      ! and face neighbour 2 is not 0
      if(any(topo%global_indices == face_nb1) .and. (face_nb2 .ne. 0)) then
         call compute_connectivity_add_connection(face_nb1, face_nb2, topo, tmp_int2d)
      endif

      ! If face neighbour 2 is local to the current rank
      ! and face neighbour 1 is not 0
      if(any(topo%global_indices == face_nb2) .and. (face_nb1 .ne. 0)) then
         call compute_connectivity_add_connection(face_nb2, face_nb1, topo, tmp_int2d)
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

    call flatten_connectivity(tmp_int2d, topo)

    call dprint ("Number of halo cells after partitioning: "//str(topo%halo_num_cells))

    topo%total_num_cells = topo%local_num_cells + topo%halo_num_cells

    call dprint ("Total number of cells (local + halo) after partitioning: "//str(topo%total_num_cells))

    !! ! First the indices of the local cells
    !! k = 1
    !! do i = 1, topo%global_num_cells
    !!   if(topo%global_partition(i) == irank) then
    !!     topo%global_indices(k) = i
    !!     call dprint("Index of local cell: topo%global_indices("//str(k)//") = "//str(i))
    !!     k = k + 1
    !!   end if 
    !! end do
    !! ! Then the indices of the halo cells
    !! k = topo%local_num_cells + 1
    !! do i = 1, size(topo%adjncy)
    !!   if (topo%global_partition(topo%adjncy(i)) /= irank) then
    !!     call dprint("Index of halo cell: topo%global_indices("//str(k)//") = "//str(i))
    !!     topo%global_indices(k) = topo%adjncy(i)
    !!     topo%global_indices(k) = int(topo%adjncy(i)) ! XXX: OK for small meshes...
    !!     k = k + 1
    !!   end if
    !! end do

    !! ! Now compute the face indices for the local vector
    !! allocate(topo%face_indices(topo%max_faces, topo%local_num_cells))
    !! do i = 1, topo%local_num_cells
    !!   k = topo%global_indices(i)
    !!   do j = 1, topo%max_faces
    !!     topo%face_indices(j,i) = topo%global_face_indices(j,k) ! XXX: currently using global faces numbers
    !!   end do
    !! end do

    !! do i = 1, size(topo%global_indices) 
    !!   call dprint("Global index "//str(i)//": "//str(int(topo%global_indices(i))))
    !! end do

    !! deallocate(topo%global_face_indices) ! No longer needed now we have the local face indices

    !! topo%num_faces = 0
    !! previous = 0

    !! ! Next, compute the number of local faces, which  
    !! ! equals the maximum value in the face_indices array
    !! topo%num_faces = maxval(topo%face_indices)

    !! call dprint("Number of local faces "//str(topo%num_faces))

    !! if(irank == 0) print*,"Global boundaries: ", topo%global_boundaries(1:10)

    ! ! Finally, allocate and compute the neighbour indices
    ! if (allocated(topo%nb_indices)) then
    !   deallocate(topo%nb_indices)
    ! end if
    ! allocate(topo%nb_indices(topo%max_faces, topo%local_num_cells))

    ! do i = 1, size(topo%xadj) - 1
    !   n = int(topo%xadj(i+1) - topo%xadj(i)) ! XXX: OK for small meshes!
    !   do j = 1, n
    !     topo%nb_indices(j, i) = int(topo%adjncy(j+topo%xadj(i)-1)) ! XXX: OK for small meshes!
    !   end do 
    !   do k = n + 1, topo%max_faces
    !     topo%nb_indices(k, i) = 0 ! Set boundaries to 0 - will be updated with correct BCs
    !   end do
    ! end do

    ! if (irank == 0) print*, "Neighbour indices ", topo%nb_indices

  end subroutine compute_connectivity

  subroutine compute_connectivity_get_local_cells(par_env, topo)

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(topology), target, intent(inout) :: topo                           !< The topology for which to compute the parition

    integer :: i
    integer :: ctr

    ! Allocate and then compute global indices
    if(allocated(topo%global_indices)) then
      deallocate(topo%global_indices)
    end if
    allocate(topo%global_indices(topo%local_num_cells))
    topo%global_indices(:) = -1 ! This will allow us to check later

    ctr = 1
    associate( irank => par_env%proc_id, &
               partition => topo%global_partition )
      do i = 1, topo%global_num_cells
        if (partition(i) == irank) then
          topo%global_indices(ctr) = i
          ctr = ctr + 1
        end if
      end do
    end associate

    if (ctr /= (topo%local_num_cells + 1)) then
      print *, "ERROR: didn't find all my cells!"
      stop
    end if

    if (minval(topo%global_indices) < 1) then
      print *, "ERROR: didn't register all cells properly!"
      stop
    end if

    if (maxval(topo%global_indices) > topo%global_num_cells) then
      print *, "ERROR: global index exceeds range!"
      stop
    end if

  end subroutine

  subroutine compute_connectivity_add_connection(face_nb1, face_nb2, topo, tmp_int2d)

    integer(ccs_int), intent(in) :: face_nb1                      !< Local cell global index
    integer(ccs_int), intent(in) :: face_nb2                      !< Neighbouring cell global index
    type(topology), target, intent(inout) :: topo        !< The topology for which to compute the partition
    integer, dimension(:, :), intent(inout) :: tmp_int2d !< Temporary connectivity array

    integer, dimension(1) :: local_index
    integer :: fctr

    local_index = findloc(topo%global_indices, face_nb1)
    if (local_index(1) <= 0) then
      print *, "ERROR: failed to find face neighbour in global indices, findloc: "
      print *, "- ANY: ", any(topo%global_indices == face_nb1)
      print *, "- local_index: ", local_index
      print *, "- Face neighbour index: ", face_nb1
      print *, "- Global indices: ", topo%global_indices
      stop
    end if
         
    fctr = tmp_int2d(local_index(1), topo%max_faces + 1) + 1 ! Increment number of faces for this cell
    tmp_int2d(local_index(1), fctr) = face_nb2               ! Store global index of neighbour cell
    tmp_int2d(local_index(1), topo%max_faces + 1) = fctr     ! Store number of faces for this cell
    topo%num_nb(local_index(1)) = fctr

  end subroutine

  !v Take the 2D connectivity graph and convert to 1D
  !  Note that cell neighbours are still globally numbered at this point.
  subroutine flatten_connectivity(tmp_int2d, topo)
    
    integer, dimension(:, :), intent(in) :: tmp_int2d
    type(topology), target, intent(inout) :: topo        !< The topology for which to compute the partition

    integer :: i, j
    integer, dimension(:), allocatable :: tmp1
    integer, dimension(:), allocatable :: tmp2

    integer :: ctr
    integer, dimension(1) :: local_idx

    topo%halo_num_cells = 0

    ctr = 1

    allocate(tmp1(topo%local_num_cells))
    tmp1(:) = -1

    do i = 1, topo%local_num_cells
      topo%xadj(i) = ctr

      ! Loop over connections of cell i
      do j = 1, tmp_int2d(i, topo%max_faces + 1)
        associate( nbidx => tmp_int2d(i, j) )
          if (.not. any(topo%global_indices == nbidx)) then
            ! Halo cell
            if (.not. any(tmp1 == nbidx)) then
              ! New halo cell
              ! Copy and extend size of halo cells buffer

              topo%halo_num_cells = topo%halo_num_cells + 1
              if (topo%halo_num_cells > size(tmp1)) then
                allocate(tmp2(size(tmp1) + topo%local_num_cells))

                tmp2(:) = -1
                tmp2(1:size(tmp1)) = tmp1(1:size(tmp1))

                deallocate(tmp1)
                allocate(tmp1, source=tmp2)
                deallocate(tmp2)
              end if

              tmp1(topo%halo_num_cells) = nbidx
            end if

            local_idx = findloc(tmp1, nbidx)
            !topo%adjncy(ctr) = local_idx(1)
          end if

          topo%adjncy(ctr) = nbidx
        end associate
            
        ctr = ctr + 1
      end do ! End j

    end do
    topo%xadj(topo%local_num_cells + 1) = ctr

    allocate(tmp2(topo%local_num_cells + topo%halo_num_cells))
    do i = 1, topo%local_num_cells
      tmp2(i) = topo%global_indices(i)
    end do
    do i = 1, topo%halo_num_cells
      tmp2(topo%local_num_cells + i) = tmp1(i)
    end do
    deallocate(topo%global_indices)
    allocate(topo%global_indices, source=tmp2)

    deallocate(tmp1)
    deallocate(tmp2)

  end subroutine

end submodule
