submodule(reordering) reordering_rcm
#include "ccs_macros.inc"

  use utils, only: debug_print
  use types, only: cell_locator, neighbour_locator

  use rcm_mod
  use rcm_kinds
  
  implicit none

contains

  !v Determine how the mesh should be reordered using bundled RCM reordering
  module subroutine get_reordering(new_indices)

    use rcm_mod
    use meshing, only: get_local_num_cells, create_cell_locator, count_neighbours, &
                       get_local_index, create_neighbour_locator, get_local_status
    use mesh_utils, only: build_adjacency_matrix

    integer(ccs_int), dimension(:), allocatable, intent(out) :: new_indices !< new indices in "to(from)" format

    integer(rcm_int), allocatable, dimension(:) :: perm, perm_inv
    integer(rcm_int) :: node_num, adj_num
    integer(rcm_int), allocatable, dimension(:) :: xadj, adjncy

    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: i

    call dprint("Reordering with RCM.")
    ! First build adjacency matrix for local cells
    call build_adjacency_matrix(xadj, adjncy)

    node_num = -1
    do i = 1, size(adjncy)
      if (adjncy(i) .gt. node_num) then
        node_num = adjncy(i)
      end if
    end do

    node_num = size(xadj) - 1
    adj_num = size(adjncy)

    allocate (perm(node_num))
    allocate (perm_inv(node_num))

    call genrcm(node_num, adj_num, xadj, adjncy, perm)
    call perm_inverse3(node_num, perm, perm_inv)

    ! Fill local indices in original ordering -> destination, i.e. to(i) => new index of cell i.
    call get_local_num_cells(local_num_cells)
    allocate (new_indices(local_num_cells))

    if (local_num_cells >= 1) then
      do i = 1, local_num_cells
        new_indices(perm(i)) = i
      end do
    end if

  end subroutine get_reordering

end submodule
