submodule(reordering) reordering_common
#include "ccs_macros.inc"

  use utils, only: exit_print, str, debug_print
  use kinds, only: ccs_int, ccs_real, ccs_err
  use types, only: cell_locator, neighbour_locator
  use meshing, only: get_local_num_cells, set_cell_location, count_neighbours, &
                     get_local_index, set_neighbour_location, get_local_status, &
                     get_centre, set_centre, &
                     get_global_index

  implicit none

  logical :: write_csr = .false.
  logical :: perform_reordering = .true.

contains

  module subroutine disable_reordering()
    perform_reordering = .false.
  end subroutine disable_reordering

  !v Cell reordering.
  !
  !  Performs a reordering of local cells and reassigns their global indices based on this new
  !  ordering.
  module subroutine reorder_cells(mesh)

    type(ccs_mesh), intent(inout) :: mesh !< the mesh to be reordered

    integer(ccs_int) :: i, wrunit
    integer(ccs_int) :: local_num_cells

    integer(ccs_int), dimension(:), allocatable :: new_indices

    if (.not. perform_reordering) then
      return
    end if

    call get_local_num_cells(mesh, local_num_cells)

    if (write_csr) then
      open (newunit=wrunit, FILE="csr_orig.txt", FORM="FORMATTED")
      do i = 1, mesh%topo%local_num_cells
        write (wrunit, *) mesh%topo%nb_indices(:, i)
      end do
      close (wrunit)
    end if

    call dprint("*********BEGIN REORDERING*****************")
    call get_reordering(mesh, new_indices)
    call dprint("---------APPLY REORDERING-----------------")
    call apply_reordering(new_indices, mesh)
    deallocate (new_indices)

    ! System indices of local indices should be contiguous
    call dprint("---------CONTIGUOUS ORDERING--------------")
    do i = 2, local_num_cells
      if (mesh%topo%global_indices(i) /= (mesh%topo%global_indices(i - 1) + 1)) then
        call error_abort("ERROR: failed global index check at local index " // str(i))
      end if
    end do
    call dprint("*********END   REORDERING*****************")

    if (write_csr) then
      open (newunit=wrunit, FILE="csr_new.txt", FORM="FORMATTED")
      do i = 1, mesh%topo%local_num_cells
        write (wrunit, *) mesh%topo%nb_indices(:, i)
      end do
      close (wrunit)
    end if

  end subroutine reorder_cells

  module subroutine apply_reordering(new_indices, mesh)

    integer(ccs_int), dimension(:), intent(in) :: new_indices !< new indices in "to(from)" format
    type(ccs_mesh), intent(inout) :: mesh                     !< the mesh to be reordered

    call dprint("         NATURAL INDICES")
    call reorder_natural_indices(new_indices, mesh)
    call dprint("         CELL CENTRES")
    call reorder_cell_centres(new_indices, mesh)
    call dprint("         CELL NEIGHBOURS")
    call reorder_neighbours(new_indices, mesh)
    
  end subroutine apply_reordering

  subroutine set_global_indices(mesh)

    use mpi

    type(ccs_mesh), intent(inout) :: mesh                     !< the mesh to be reordered

    integer(ccs_int) :: i
    integer(ccs_int) :: idxg
    integer(ccs_int) :: offset
    
    integer(ccs_int), dimension(:), allocatable :: global_indices

    integer(ccs_int) :: nrank
    integer(ccs_int) :: rank_idx
    integer(ccs_err) :: ierr
    
    if (allocated(mesh%topo%global_indices)) then
      deallocate(mesh%topo%global_indices)
    end if
    allocate(mesh%topo%global_indices(mesh%topo%total_num_cells))

    call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
    rank_idx = nrank + 1 ! For indexing into stuff...
    offset = int(mesh%topo%vtxdist(rank_idx), ccs_int) ! Everything else is ccs_int...
    if (mesh%topo%vtxdist(1) == 0) then
      ! Using C numbering
      offset = offset + 1
    end if
    do i = 1, mesh%topo%local_num_cells
      mesh%topo%global_indices(i) = i + (offset - 1)
    end do
    
    allocate(global_indices(mesh%topo%global_num_cells))

    global_indices(:) = 0

    do i = 1, mesh%topo%local_num_cells
      idxg = mesh%topo%natural_indices(i)
      global_indices(idxg) = mesh%topo%global_indices(i)
    end do
    call MPI_Allreduce(MPI_IN_PLACE, global_indices, mesh%topo%global_num_cells, &
                       MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    do i = mesh%topo%local_num_cells + 1, mesh%topo%total_num_cells
      idxg = mesh%topo%natural_indices(i)
      mesh%topo%global_indices(i) = global_indices(idxg)
    end do
    
    deallocate(global_indices)
    
  end subroutine set_global_indices

  !v Store the natural indices of the problem in the new ordering.
  !
  !  Note that up to this point the "natural" indices are stored in the global indices - these will
  !  be replaced by the global indexing of the linear system.
  !
  !  Halo cells are left in place, therefore only the natural indices of local cells need to be
  !  reordered. Note, however, that the global (linear system) index of halo cells does need to be 
  !  updated by the call to set_global_indices.
  subroutine reorder_natural_indices(new_indices, mesh)

    integer(ccs_int), dimension(:), intent(in) :: new_indices !< new indices in "to(from)" format
    type(ccs_mesh), intent(inout) :: mesh                     !< the mesh to be reordered

    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: i
    integer(ccs_int) :: idx_new
    integer(ccs_int) :: idxg
    type(cell_locator) :: loc_p
    
    ! Only need to order local cells
    call get_local_num_cells(mesh, local_num_cells)

    associate(total_num_cells => mesh%topo%total_num_cells)
      allocate(mesh%topo%natural_indices(total_num_cells))
    end associate

    ! Copy the global indices into natural indices
    mesh%topo%natural_indices(:) = mesh%topo%global_indices(:)

    ! Apply reordering on the local natural indices
    do i = 1, local_num_cells
      call set_cell_location(mesh, i, loc_p)
      call get_global_index(loc_p, idxg)
      
      idx_new = new_indices(i)
      mesh%topo%natural_indices(idx_new) = idxg
    end do

    ! Set global indices to linear system
    call set_global_indices(mesh)
    
  end subroutine reorder_natural_indices
  
  subroutine reorder_cell_centres(new_indices, mesh)

    integer(ccs_int), dimension(:), intent(in) :: new_indices !< new indices in "to(from)" format
    type(ccs_mesh), intent(inout) :: mesh                     !< the mesh to be reordered

    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: i
    integer(ccs_int) :: idx_new
    type(cell_locator) :: loc_p
    
    real(ccs_real), dimension(:, :), allocatable :: x

    call get_local_num_cells(mesh, local_num_cells)

    allocate (x(3, local_num_cells))

    do i = 1, local_num_cells
      idx_new = new_indices(i)
      call set_cell_location(mesh, i, loc_p)
      call get_centre(loc_p, x(:, idx_new))
    end do
    do i = 1, local_num_cells
      call set_cell_location(mesh, i, loc_p)
      call set_centre(loc_p, x(:, i))
    end do

    deallocate (x)

  end subroutine reorder_cell_centres

  subroutine reorder_neighbours(new_indices, mesh)

    integer(ccs_int), dimension(:), intent(in) :: new_indices !< new indices in "to(from)" format
    type(ccs_mesh), intent(inout) :: mesh                     !< the mesh to be reordered

    integer(ccs_int) :: local_num_cells
    
    integer(ccs_int) :: i, j
    integer(ccs_int) :: idx_tmp
    integer(ccs_int) :: idx_new
    
    integer(ccs_int), dimension(:, :), allocatable :: idx_nb
    integer(ccs_int), dimension(:), allocatable :: num_nb

    type(cell_locator) :: loc_p
    integer(ccs_int) :: nnb
    
    call get_local_num_cells(mesh, local_num_cells)

    allocate (idx_nb, mold=mesh%topo%nb_indices)

    idx_nb(:, :) = mesh%topo%nb_indices(:, :)
    do i = 1, local_num_cells
      call set_cell_location(mesh, i, loc_p)
      call count_neighbours(loc_p, nnb)

      ! Get new /local/ index of neighbours, note only local cells are reordered, the local
      ! index of halo cells remains unchanged.
      do j = 1, nnb
        idx_tmp = mesh%topo%nb_indices(j, i)
        if ((idx_tmp > 0) .and. (idx_tmp <= local_num_cells)) then
          idx_new = new_indices(idx_tmp)
          idx_tmp = new_indices(i)
          idx_nb(j, idx_tmp) = idx_new
        else
          idx_tmp = new_indices(i)
          idx_nb(j, idx_tmp) = mesh%topo%nb_indices(j, i)
        end if
      end do
    end do
    do i = 1, local_num_cells
      mesh%topo%nb_indices(:, i) = idx_nb(:, i)
    end do

    deallocate (idx_nb)

    allocate(num_nb(local_num_cells))
    num_nb(new_indices(:)) = mesh%topo%num_nb(:)
    mesh%topo%num_nb(:) = num_nb(:)
    deallocate(num_nb)

  end subroutine reorder_neighbours
  
  module subroutine bandwidth(mesh)

    type(ccs_mesh), intent(in) :: mesh !< the mesh to evaluate

    integer(ccs_int) :: bw, bw_max, bw_rowmax
    real(ccs_real) :: bw_avg
    integer(ccs_int) :: idx_p, idx_nb

    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb

    logical :: is_local

    integer(ccs_int) :: i, j, nnb
    integer(ccs_int) :: local_num_cells

    bw_max = 0
    bw_avg = 0.0

    call get_local_num_cells(mesh, local_num_cells)
    do i = 1, local_num_cells
      call set_cell_location(mesh, i, loc_p)
      call count_neighbours(loc_p, nnb)
      call get_local_index(loc_p, idx_p)
      bw_rowmax = 0
      do j = 1, nnb
        call set_neighbour_location(loc_p, j, loc_nb)
        call get_local_status(loc_nb, is_local)
        if (is_local) then
          call get_local_index(loc_nb, idx_nb)
          bw = abs(idx_nb - idx_p)
          bw_max = max(bw_max, bw)
          bw_rowmax = max(bw_rowmax, bw)
        end if
      end do
      bw_avg = bw_avg + bw_rowmax
    end do
    if (local_num_cells > 0) then
      bw_avg = bw_avg / local_num_cells
      print *, "Bandwidth: ", bw_max, bw_avg
    end if

  end subroutine bandwidth

end submodule
