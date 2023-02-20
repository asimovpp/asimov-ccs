submodule(reordering) reordering_rcm
  
  use types, only: cell_locator, neighbour_locator

  implicit none

contains

  module subroutine get_reordering(mesh, new_indices)

    use rcm_mod
    use meshing, only: get_local_num_cells, set_cell_location, count_neighbours, &
                       get_local_index, set_neighbour_location, get_local_status

    type(ccs_mesh), intent(in) :: mesh
    integer(ccs_int), dimension(:), allocatable, intent(out) :: new_indices

    integer(ccs_int), allocatable, dimension(:) :: perm, perm_inv
    integer(ccs_int) :: node_num, adj_num 
    integer(ccs_int), allocatable, dimension(:) :: xadj, adjncy

    integer(ccs_int) :: local_num_cells

    integer(ccs_int) :: i, j, nnb
    integer(ccs_int) :: ctr
    integer(ccs_int) :: idx 
    logical :: cell_local
    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb

    integer(ccs_int) :: idx_new
    
    ! First build adjacency matrix for local cells
    call get_local_num_cells(mesh, local_num_cells)

    allocate(xadj(0))
    allocate(adjncy(0))
    ctr = 1
    xadj = [xadj, ctr]
    do i = 1, local_num_cells
       call set_cell_location(mesh, i, loc_p)
       call count_neighbours(loc_p, nnb)
       do j = 1, nnb
          call set_neighbour_location(loc_p, j, loc_nb)
          call get_local_status(loc_nb, cell_local)
          if (cell_local) then
             call get_local_index(loc_nb, idx)
             adjncy = [adjncy, idx]
             ctr = ctr + 1
          end if
       end do      
       xadj = [xadj, ctr]
    end do
    
    node_num = -1
    do i = 1, size(adjncy)
      if (adjncy(i) .gt. node_num) then
        node_num = adjncy(i)
      end if
    end do

    node_num = size(xadj) - 1
    adj_num = size(adjncy) 

    print *, "node_num", node_num, "adj_num", adj_num 

    allocate(perm(node_num))
    allocate(perm_inv(node_num))

    !call genrcm(node_num, adj_num, mesh%topo%xadj, mesh%topo%adjncy, perm)  
    ! create copies of xadj and adjncy and convert them to integers from longs
    ! because genrcm expects integers, but parhip (above) expects longs
    ! call adj_show(node_num, adj_num, xadj, adjncy)
    call genrcm(node_num, adj_num, xadj, adjncy, perm)  
    call perm_inverse3(node_num, perm, perm_inv)
    ! call i4vec_print ( node_num, perm, '  The RCM permutation:' )
    ! call adj_perm_show(node_num, adj_num, xadj, adjncy, perm, perm_inv )

    
    ! Fill local indices in original ordering -> destination, i.e. to(i) => new index of cell i.
    allocate(new_indices(local_num_cells))

    if (local_num_cells >= 1) then
       do i = 1, local_num_cells
         new_indices(perm(i)) = i
       end do
     end if

  end subroutine get_reordering

end submodule
