submodule(reordering) reordering_common
#include "ccs_macros.inc"

  use utils, only: exit_print, str, debug_print
  use kinds, only: ccs_real, ccs_err
  use types, only: cell_locator, neighbour_locator
  use meshing, only: get_local_num_cells, set_cell_location, count_neighbours, &
                     get_local_index, set_neighbour_location, get_local_status, &
                     get_centre, set_centre

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
  !  ordering - assumes a contiguous numbering of the processor's partition, i.e. proc N has global
  !  indices g0 - gN. Following reordering an update is required to inform other processors about the
  !  new global indices of their halo cells.
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
    call apply_reordering(new_indices, mesh)
    deallocate (new_indices)

    ! Global indices of local indices should be contiguous
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

    use mpi

    integer(ccs_int), dimension(:), intent(in) :: new_indices !< new indices in "to(from)" format
    type(ccs_mesh), intent(inout) :: mesh                     !< the mesh to be reordered

    integer(ccs_err) :: ierr

    integer(ccs_int) :: local_num_cells

    integer(ccs_int) :: i

    integer(ccs_int), dimension(:), allocatable :: new_global_ordering

    integer(ccs_int) :: start_global
    integer(ccs_int) :: idxg, idx_new

    call get_local_num_cells(mesh, local_num_cells)

    !! Get global reordering
    allocate (new_global_ordering(mesh%topo%global_num_cells))
    new_global_ordering(:) = 0
    start_global = -1
    if (local_num_cells >= 1) then
      start_global = mesh%topo%global_indices(1)

      do i = 1, local_num_cells
        idx_new = new_indices(i)
        idxg = i + (start_global - 1)
        new_global_ordering(idxg) = idx_new + (start_global - 1)
      end do
    end if

    call MPI_Allreduce(MPI_IN_PLACE, new_global_ordering, mesh%topo%global_num_cells, &
                       MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    !! Apply reordering

    ! Set global indexing
    do i = 1, local_num_cells
      mesh%topo%global_indices(i) = i + (start_global - 1)
    end do
    do i = local_num_cells + 1, mesh%topo%total_num_cells
      idxg = mesh%topo%global_indices(i)
      mesh%topo%global_indices(i) = new_global_ordering(idxg)
    end do
    deallocate (new_global_ordering)

    call reorder_cell_centres(new_indices, mesh)
    call reorder_neighbours(new_indices, mesh)
    
  end subroutine apply_reordering

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
    bw_avg = bw_avg / local_num_cells
    print *, "Bandwidth: ", bw_max, bw_avg

  end subroutine bandwidth

end submodule
