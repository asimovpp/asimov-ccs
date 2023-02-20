submodule(reordering) reordering_common
  
  use kinds, only: ccs_real, ccs_err
  use types, only: cell_locator, neighbour_locator

  implicit none

contains
  !v Cell reordering.
  !
  !  Performs a reordering of local cells and reassigns their global indices based on this new
  !  ordering - assumes a contiguous numbering of the processor's partition, i.e. proc N has global
  !  indices g0 - gN. Following reordering an update is required to inform other processors about the
  !  new global indices of their halo cells.
  module subroutine reorder_cells(mesh)
    
    use meshing, only: get_local_num_cells

    type(ccs_mesh), intent(inout) :: mesh

    integer(ccs_int) :: i
    integer(ccs_int) :: local_num_cells
    
    integer(ccs_int), dimension(:), allocatable :: new_indices

    call get_local_num_cells(mesh, local_num_cells)
    
    open(unit=2031, FILE="csr_orig.txt", FORM="FORMATTED")
    do i = 1, mesh%topo%local_num_cells
       write(2031, *) mesh%topo%nb_indices(:, i)
    end do
    close(2031)

    ! print *, "**************************"
    ! call get_reordering(mesh, new_indices)
    ! print *, new_indices
    print *, "**************************"
    call get_reordering(mesh, new_indices)
    ! print *, new_indices
    print *, "**************************"
    call apply_reordering(new_indices, mesh)
    deallocate(new_indices)

    ! Global indices of local indices should be contiguous
    do i = 2, local_num_cells
       if (mesh%topo%global_indices(i) /= (mesh%topo%global_indices(i - 1) + 1)) then
          print *, "ERROR: failed global index check at local index ", i
          stop
       end if
    end do
    
    open(unit=2032, FILE="csr_new.txt", FORM="FORMATTED")
    do i = 1, mesh%topo%local_num_cells
       write(2032, *) mesh%topo%nb_indices(:, i)
    end do
    close(2032)
    
  end subroutine reorder_cells
  ! Given a reordering, apply it.
  module subroutine apply_reordering(new_indices, mesh)

    use mpi
    use meshing, only: get_local_num_cells, set_cell_location, count_neighbours, &
                       get_local_index, set_neighbour_location, get_local_status, &
                       get_centre, set_centre
    
    integer(ccs_int), dimension(:), intent(in) :: new_indices
    type(ccs_mesh), intent(inout) :: mesh

    integer(ccs_err) :: ierr
    
    integer(ccs_int) :: local_num_cells

    type(cell_locator) :: loc_p
    integer(ccs_int) :: i, j

    integer(ccs_int), dimension(:), allocatable :: new_global_ordering

    integer(ccs_int) :: start_global
    integer(ccs_int) :: idxg, idx_new

    real(ccs_real), dimension(:, :), allocatable :: x
    integer(ccs_int), dimension(:, :), allocatable :: idx_nb

    call get_local_num_cells(mesh, local_num_cells)

    !! Get global reordering
    allocate(new_global_ordering(mesh%topo%global_num_cells))
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

     call MPI_Allreduce(MPI_IN_PLACE, new_global_ordering, mesh%topo%global_num_cells, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

     !! Apply reordering

     ! Set global indexing
     do i = 1, local_num_cells
       mesh%topo%global_indices(i) = i + (start_global - 1)
     end do
     do i = local_num_cells + 1, mesh%topo%total_num_cells
       idxg = mesh%topo%global_indices(i)
       mesh%topo%global_indices(i) = new_global_ordering(idxg)
     end do
     deallocate(new_global_ordering)
    
     ! Reorder cell centres
     allocate(x(3, local_num_cells))
     do i = 1, local_num_cells
       idx_new = new_indices(i)
       call set_cell_location(mesh, i, loc_p)
       call get_centre(loc_p, x(:, idx_new))
     end do
     do i = 1, local_num_cells
       call set_cell_location(mesh, i, loc_p)
       call set_centre(loc_p, x(:, i))
     end do
     deallocate(x)

     ! Reorder neighbours
     allocate(idx_nb(4, local_num_cells))
     idx_nb(:,:) = mesh%topo%nb_indices(:,:)
     do i = 1, local_num_cells ! First update the neighbour copy
       do j = 1, 4
         idxg = mesh%topo%nb_indices(j, i)
         if ((idxg > 0) .and. (idxg <= local_num_cells)) then
           idx_new = new_indices(idxg)
           idxg = new_indices(i)
           idx_nb(j, idxg) = idx_new
         else
           idxg = new_indices(i)
           idx_nb(j, idxg) = mesh%topo%nb_indices(j, i)
         end if
       end do
     end do
     do i = 1, local_num_cells
       mesh%topo%nb_indices(:, i) = idx_nb(:, i)
     end do

    deallocate(idx_nb)
     
  end subroutine apply_reordering
  
  module subroutine bandwidth(mesh)
    use meshing, only: get_local_num_cells, set_cell_location, count_neighbours, &
                       get_local_index, set_neighbour_location, get_local_status
                       

    type(ccs_mesh), intent(in) :: mesh

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
