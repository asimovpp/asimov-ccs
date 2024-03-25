submodule(reordering) reordering_common
#include "ccs_macros.inc"

  use utils, only: exit_print, str, debug_print
  use kinds, only: ccs_int, ccs_real, ccs_err
  use types, only: cell_locator, neighbour_locator
  use parallel_types, only: parallel_environment
  use meshing, only: get_local_num_cells, create_cell_locator, count_neighbours, &
                     get_local_index, create_neighbour_locator, &
                     get_boundary_status, get_local_status, &
                     get_centre, set_centre, &
                     get_global_index, set_global_index, &
                     get_natural_index, &
                     get_global_num_cells, &
                     set_natural_index, &
                     get_total_num_cells, &
                     get_vert_per_cell
  use parallel, only: create_shared_array, destroy_shared_array, is_root

  implicit none

  logical :: write_csr = .false.

contains

  !v Cell reordering.
  !
  !  Performs a reordering of local cells and reassigns their global indices based on this new
  !  ordering.
  module subroutine reorder_cells(par_env, shared_env, mesh)

    class(parallel_environment), intent(in) :: par_env !< The parallel environment
    class(parallel_environment), intent(in) :: shared_env !< The parallel environment
    type(ccs_mesh), intent(inout) :: mesh              !< the mesh to be reordered

    integer(ccs_int) :: i, wrunit
    integer(ccs_int) :: local_num_cells

    integer(ccs_int), dimension(:), allocatable :: new_indices

    call get_local_num_cells(local_num_cells)

    if (write_csr .and. is_root(par_env)) then
      open (newunit=wrunit, FILE="csr_orig.txt", FORM="FORMATTED")
      do i = 1, mesh%topo%local_num_cells
        write (wrunit, *) mesh%topo%nb_indices(:, i)
      end do
      close (wrunit)
    end if

    call dprint("*********BEGIN REORDERING*****************")
    call get_reordering(new_indices)
    call dprint("---------APPLY REORDERING-----------------")
    call apply_reordering(new_indices, par_env, shared_env, mesh)
    deallocate (new_indices)

    ! System indices of local indices should be contiguous
    call dprint("---------CONTIGUOUS ORDERING--------------")
    do i = 2, local_num_cells
      if (mesh%topo%global_indices(i) /= (mesh%topo%global_indices(i - 1) + 1)) then
        call error_abort("ERROR: failed global index check at local index " // str(i))
      end if
    end do
    call dprint("*********END   REORDERING*****************")

    if (write_csr .and. is_root(par_env)) then
      open (newunit=wrunit, FILE="csr_new.txt", FORM="FORMATTED")
      do i = 1, mesh%topo%local_num_cells
        write (wrunit, *) mesh%topo%nb_indices(:, i)
      end do
      close (wrunit)
    end if

  end subroutine reorder_cells

  module subroutine apply_reordering(new_indices, par_env, shared_env, mesh)

    integer(ccs_int), dimension(:), intent(in) :: new_indices !< new indices in "to(from)" format
    class(parallel_environment), intent(in) :: par_env        !< The parallel environment
    class(parallel_environment), intent(in) :: shared_env        !< The parallel environment
    type(ccs_mesh), intent(inout) :: mesh                     !< the mesh to be reordered

    call dprint("         NATURAL INDICES")
    call reorder_natural_indices(new_indices, par_env, mesh)
    call dprint("         CELL NEIGHBOURS")
    call reorder_neighbours(new_indices, mesh)
    call dprint("         CELL FACES")
    call reorder_faces(new_indices, mesh)
    call dprint("         CELL VERTICES")
    call reorder_vertices(shared_env, mesh)

  end subroutine apply_reordering

  subroutine set_global_indices(par_env, mesh)

    use mpi

    use parallel_types_mpi, only: parallel_environment_mpi

    class(parallel_environment), intent(in) :: par_env        !< The parallel environment
    type(ccs_mesh), intent(inout) :: mesh                     !< the mesh to be reordered

    integer(ccs_int) :: i
    integer(ccs_int) :: idxg
    integer(ccs_int) :: idxn
    integer(ccs_int) :: offset

    integer(ccs_int), dimension(:), allocatable :: global_indices

    integer(ccs_err) :: ierr

    integer(ccs_int) :: global_num_cells
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: total_num_cells
    type(cell_locator) :: loc_p

    call get_total_num_cells(total_num_cells)
    if (allocated(mesh%topo%global_indices)) then
      deallocate (mesh%topo%global_indices)
    end if
    allocate (mesh%topo%global_indices(total_num_cells))

    ! Compute the global offset, this should start from 1.
    call get_global_offset(par_env, offset)

    ! Apply local (contiguous) global numbering
    call get_local_num_cells(local_num_cells)
    do i = 1, local_num_cells
      call create_cell_locator(i, loc_p)
      call set_global_index(i + (offset - 1), loc_p)
    end do

    ! Determine the new global index of halo cells.
    ! The easiest way to do this is a global array with the new global indices in original ordering, i.e. to(from).
    call get_global_num_cells(global_num_cells)
    allocate (global_indices(global_num_cells))

    global_indices(:) = 0

    do i = 1, local_num_cells
      call create_cell_locator(i, loc_p)
      call get_natural_index(loc_p, idxn) ! where the cell was in the original ordering
      call get_global_index(loc_p, idxg)
      global_indices(idxn) = idxg
    end do

    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_Allreduce(MPI_IN_PLACE, global_indices, mesh%topo%global_num_cells, &
                         MPI_INTEGER, MPI_SUM, par_env%comm, ierr)
    class default
      call error_abort("Unsupported parallel environment!")
    end select

    do i = local_num_cells + 1, total_num_cells
      call create_cell_locator(i, loc_p)
      call get_natural_index(loc_p, idxn)
      idxg = global_indices(idxn)
      call set_global_index(idxg, loc_p)
    end do

    deallocate (global_indices)

  end subroutine set_global_indices

  !v Store the natural indices of the problem in the new ordering.
  !
  !  Note that up to this point the "natural" indices are stored in the global indices - these will
  !  be replaced by the global indexing of the linear system.
  !
  !  Halo cells are left in place, therefore only the natural indices of local cells need to be
  !  reordered. Note, however, that the global (linear system) index of halo cells does need to be
  !  updated by the call to set_global_indices.
  subroutine reorder_natural_indices(new_indices, par_env, mesh)

    integer(ccs_int), dimension(:), intent(in) :: new_indices !< new indices in "to(from)" format
    class(parallel_environment), intent(in) :: par_env        !< The parallel environment
    type(ccs_mesh), intent(inout) :: mesh                     !< the mesh to be reordered

    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: total_num_cells
    integer(ccs_int) :: i
    integer(ccs_int) :: idx_new
    integer(ccs_int) :: idxg
    type(cell_locator) :: loc_p

    ! Only need to order local cells
    call get_local_num_cells(local_num_cells)
    call get_total_num_cells(total_num_cells)
    allocate (mesh%topo%natural_indices(total_num_cells))
    mesh%topo%natural_indices(:) = 0 ! For checking

    ! Apply reordering on the local natural indices
    do i = 1, local_num_cells
      call create_cell_locator(i, loc_p)
      call get_global_index(loc_p, idxg)

      idx_new = new_indices(i)
      call create_cell_locator(idx_new, loc_p)
      call set_natural_index(idxg, loc_p)
    end do

    ! Copy halo global indices -> natural indices
    do i = local_num_cells + 1, total_num_cells
      call create_cell_locator(i, loc_p)
      call get_global_index(loc_p, idxg)
      call set_natural_index(idxg, loc_p)
    end do

    ! Set global indices to linear system
    call set_global_indices(par_env, mesh)

  end subroutine reorder_natural_indices

  subroutine reorder_neighbours(new_indices, mesh)

    integer(ccs_int), dimension(:), intent(in) :: new_indices !< new indices in "to(from)" format
    type(ccs_mesh), intent(inout) :: mesh                     !< the mesh to be reordered

    integer(ccs_int) :: local_num_cells

    integer(ccs_int) :: i, j
    integer(ccs_int) :: idx_tmp
    integer(ccs_int) :: idx_new

    integer(ccs_int), dimension(:, :), allocatable :: idx_nb

    type(cell_locator) :: loc_p
    integer(ccs_int) :: nnb

    type(neighbour_locator) :: loc_nb
    logical :: is_boundary
    logical :: is_local

    call get_local_num_cells(local_num_cells)

    allocate (idx_nb, source=mesh%topo%nb_indices)

    do i = 1, local_num_cells
      call create_cell_locator(i, loc_p)
      call count_neighbours(loc_p, nnb)

      ! Get new /local/ index of neighbours, note only local cells are reordered, the local
      ! index of halo cells remains unchanged.
      do j = 1, nnb
        call create_neighbour_locator(loc_p, j, loc_nb)
        call get_boundary_status(loc_nb, is_boundary)
        call get_local_status(loc_nb, is_local)
        if ((.not. is_boundary) .and. is_local) then
          call get_local_index(loc_nb, idx_tmp)
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

    ! Reorder neighbour counts
    mesh%topo%num_nb(new_indices(:)) = mesh%topo%num_nb(:)

  end subroutine reorder_neighbours

  subroutine reorder_faces(new_indices, mesh)

    integer(ccs_int), dimension(:), intent(in) :: new_indices !< new indices in "to(from)" format
    type(ccs_mesh), intent(inout) :: mesh                     !< the mesh to be reordered

    mesh%topo%face_indices(:, new_indices(:)) = mesh%topo%face_indices(:, :)

  end subroutine

  !> Extract the vertex indices from the global data (restricting to local data)
  subroutine reorder_vertices(shared_env, mesh)

    class(parallel_environment), intent(in) :: shared_env
    type(ccs_mesh), intent(inout) :: mesh                     !< the mesh to be reordered

    integer(ccs_int) :: vert_per_cell
    integer(ccs_err) :: local_num_cells

    integer(ccs_int) :: i

    type(cell_locator) :: loc_p
    integer(ccs_int) :: natural_index

    call get_vert_per_cell(vert_per_cell)
    call get_local_num_cells(local_num_cells)
    allocate (mesh%topo%loc_global_vertex_indices(vert_per_cell, local_num_cells))

    ! Extract vertex indices of local cells from the global array
    do i = 1, local_num_cells
      call create_cell_locator(i, loc_p)
      call get_natural_index(loc_p, natural_index)
      mesh%topo%loc_global_vertex_indices(:, i) = mesh%topo%global_vertex_indices(:, natural_index)
    end do

    call destroy_shared_array(shared_env, mesh%topo%global_vertex_indices, mesh%topo%global_vertex_indices_window)


  end subroutine

  module subroutine print_bandwidth(par_env)

    use mpi
    
    use case_config, only: compute_bwidth
    use parallel_types_mpi, only: parallel_environment_mpi

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    
    integer(ccs_int) :: bw_max
    real(ccs_real) :: bw_avg

    integer(ccs_int) :: ierr
    real(ccs_real) :: sum_bw_avg

    if (.not. compute_bwidth) then
      return
    end if

    call compute_bandwidth(bw_max, bw_avg)
    
    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_Allreduce(MPI_IN_PLACE, bw_max, 1, MPI_INTEGER, MPI_MAX, par_env%comm, ierr)
      call MPI_Allreduce(bw_avg, sum_bw_avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, par_env%comm, ierr)
      if(is_root(par_env)) then 
        print *, "Bandwidth: ", bw_max, sum_bw_avg/par_env%num_procs
      end if

    class default
      call error_abort("Unsupported parallel environment!")
    end select
    
  end subroutine

  subroutine compute_bandwidth(bw_max, bw_avg)

    integer(ccs_int), intent(out) :: bw_max
    real(ccs_real), intent(out) :: bw_avg

    integer(ccs_int) :: bw, bw_rowmax
    integer(ccs_int) :: idx_p, idx_nb

    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb

    logical :: is_local

    integer(ccs_int) :: i, j, nnb
    integer(ccs_int) :: local_num_cells

    bw_max = 0
    bw_avg = 0.0

    call get_local_num_cells(local_num_cells)
    do i = 1, local_num_cells
      call create_cell_locator(i, loc_p)
      call count_neighbours(loc_p, nnb)
      call get_local_index(loc_p, idx_p)
      bw_rowmax = 0
      do j = 1, nnb
        call create_neighbour_locator(loc_p, j, loc_nb)
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

  end subroutine compute_bandwidth

  ! Get the cell distribution across all processors in rank order, and compute my offset.
  subroutine get_global_offset(par_env, offset)

    use mpi

    use parallel_types_mpi, only: parallel_environment_mpi

    class(parallel_environment), intent(in) :: par_env
    integer(ccs_int), intent(out) :: offset

    integer(ccs_int) :: nproc
    integer(ccs_int) :: par_idx
    integer(ccs_int), dimension(:), allocatable :: cell_counts

    integer(ccs_int) :: local_num_cells

    integer(ccs_int) :: i
    integer(ccs_err) :: ierr

    select type (par_env)
    type is (parallel_environment_mpi)
      nproc = par_env%num_procs
      par_idx = par_env%proc_id + 1 ! MPI is C-indexed
    class default
      call error_abort("Unsupported parallel environment!")
    end select

    allocate (cell_counts(nproc))

    call get_local_num_cells(local_num_cells)
    cell_counts(:) = 0
    cell_counts(par_idx) = local_num_cells

    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_Allreduce(MPI_IN_PLACE, cell_counts, nproc, &
                         MPI_INTEGER, MPI_SUM, par_env%comm, &
                         ierr)
    class default
      call error_abort("Unsupported parallel environment!")
    end select

    offset = 1
    do i = 1, par_idx - 1
      offset = offset + cell_counts(i)
    end do

    deallocate (cell_counts)

  end subroutine

end submodule
