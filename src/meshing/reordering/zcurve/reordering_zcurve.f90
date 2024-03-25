submodule(reordering) reordering_zcurve
#include "ccs_macros.inc"

  use utils, only: debug_print
  use types, only: cell_locator, neighbour_locator
  use kinds, only: ccs_long
  use parallel, only: is_root

  implicit none

contains

  !v Determine how the mesh should be reordered according to the Morton space filling curve
  module subroutine get_reordering(new_indices)

    use mortif, only: demorton2D  
    use meshing, only: get_local_num_cells, create_cell_locator, count_neighbours, &
                       get_local_index, create_neighbour_locator, get_local_status
    use mesh_utils, only: build_adjacency_matrix

    integer(ccs_int), dimension(:), allocatable, intent(out) :: new_indices !< new indices in "to(from)" format

    integer(ccs_int), allocatable, dimension(:) :: xadj, adjncy
    integer(ccs_int) :: local_num_cells

    integer(ccs_int) :: i
    integer(ccs_int) :: ctr

    integer(ccs_long) :: mcode
    integer(ccs_long) :: max_mcode
    integer(ccs_int) :: row, col
    logical :: isfound
    integer(ccs_int), dimension(:), allocatable :: perm 

    call dprint("Reordering with z-curve.")

    ! First build adjacency matrix for local cells
    call build_adjacency_matrix(xadj, adjncy)

    ! Space filling curve needs to be big enough to cover the adjacency matrix.
    ! Deal with non-power-of-2-cases
    call get_local_num_cells(local_num_cells)
    max_mcode = 1
    do while (max_mcode < local_num_cells)
      max_mcode = max_mcode * 2
    end do
    
    allocate (perm(local_num_cells))
    perm(:) = -1
    ctr = 1

    ! Now traverse space filling curve and generate new indices
    do mcode = 0, max_mcode**2 - 1
      call demorton2D(code=mcode, i=col, j=row)
      ! morton is calculated as 0 index, so I should add 1 to each dim when addressing Fortran arrays
      col = col + 1
      row = row + 1

      ! the zcurve traverses a bigger area than encompassed by the adjacency matrix, so check if we are inside the adj matrix
      if (max(col, row) <= local_num_cells) then

        ! find if the (col, row) cell is in the adjacency matrix
        isfound = .false.
        ! the col==row non-zeros are not explicitly in the adj matrix
        if (col .eq. row) then
          isfound = .true.
        else if (any(col .eq. adjncy(xadj(row):xadj(row+1)-1))) then 
          isfound = .true.
        end if

        if (isfound) then
          ! make sure each column is permuted only once
          !if (.not. any(col == perm(1:ctr))) then !TODO: check this optimisation
          if (.not. any(col == perm(:))) then
            perm(ctr) = col
            ctr = ctr + 1
          end if
        end if
        
      end if
    end do

    ! Transform indices to format: to(i) => new index of cell i.
    allocate (new_indices(local_num_cells))
    if (local_num_cells >= 1) then
      do i = 1, local_num_cells
        new_indices(perm(i)) = i
      end do
    end if

  end subroutine get_reordering

end submodule
