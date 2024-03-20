submodule(partitioning) partitioning_common
#include "ccs_macros.inc"

  use kinds, only: ccs_int, ccs_err, ccs_real
  use types, only: topology, graph_connectivity, cell_locator, neighbour_locator
  use utils, only: str, debug_print, exit_print
  use parallel_types_mpi, only: parallel_environment_mpi
  use mesh_utils, only: count_mesh_faces, set_cell_face_indices
  use meshing, only: set_local_num_cells, get_local_num_cells, get_global_num_cells, &
                     get_halo_num_cells, set_halo_num_cells, &
                     get_global_num_faces, &
                     get_total_num_cells, set_total_num_cells, &
                     set_num_faces, &
                     get_max_faces, &
                     create_cell_locator, create_neighbour_locator, &
                     get_global_index, &
                     get_local_index, set_local_index, &
                     set_mesh_object, nullify_mesh_object, &
                     set_topo_object, nullify_topo_object
  use parallel, only: is_root, is_valid, create_shared_array, destroy_shared_array, sync

  implicit none

contains

  ! Compute the new topology connectivity after partitioning
  module subroutine compute_connectivity(par_env, shared_env, mesh)

    use mpi

    class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The global parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared parallel environment
    type(ccs_mesh), target, intent(inout) :: mesh                              !< The mesh for which to compute the parition

    ! Local variables
    integer(ccs_int) :: irank ! MPI rank ID
    integer(ccs_int) :: isize ! Size of MPI world
    integer(ccs_int) :: i
    integer(ccs_int) :: partition_rank
    
    irank = par_env%proc_id
    isize = par_env%num_procs

    ! Get global indices of local cells
    call compute_connectivity_get_local_cells(par_env, mesh)

    ! Recompute vtxdist array based on the new partition
    associate(global_partition => mesh%topo%graph_conn%global_partition, &
              vtxdist => mesh%topo%graph_conn%vtxdist)
      vtxdist(:) = 0

      ! Store counts at rank + 1
      do i = 1, size(global_partition)
        partition_rank = int(global_partition(i), ccs_int) + 1 ! Ranks are C-indexed...
        vtxdist(partition_rank + 1) = vtxdist(partition_rank + 1) + 1
      end do

      ! Compute distribution offsets
      vtxdist(1) = 1
      do i= 2, isize + 1
        vtxdist(i) = vtxdist(i) + vtxdist(i - 1)
      end do
    end associate
    
    if (irank == 0) then
      do i = 1, isize + 1
        call dprint("new vtxdist(" // str(i) // "): " // str(int(mesh%topo%graph_conn%vtxdist(i))))
      end do
    end if

    call compute_face_connectivity(par_env, shared_env, mesh)

  end subroutine compute_connectivity

  subroutine compute_face_connectivity(par_env, shared_env, mesh)

    use iso_fortran_env, only: int32

    class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared parallel environment
    type(ccs_mesh), target, intent(inout) :: mesh                              !< The mesh for which to compute the parition

    ! Local variables
    integer(ccs_int), dimension(:, :), allocatable :: tmp_int2d ! Temporary 2D integer array
    integer(ccs_int), dimension(:, :), pointer :: cell_faces
    integer(ccs_int), dimension(:), allocatable :: cell_faces_counters
    integer :: cell_faces_window
    integer(ccs_int) :: irank ! MPI rank ID
    integer(ccs_int) :: i
    integer(ccs_int) :: j
    integer(ccs_int) :: face_nb1
    integer(ccs_int) :: face_nb2
    integer(ccs_int) :: num_connections
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: halo_num_cells
    integer(ccs_int) :: total_num_cells
    integer(ccs_int) :: global_num_cells
    integer(ccs_int) :: global_num_faces
    integer(ccs_int) :: max_faces
    
    integer(ccs_int) :: global_index_p
    integer(ccs_int) :: face
    integer(ccs_int) :: neighbour

    irank = par_env%proc_id

    call get_local_num_cells(local_num_cells)

    ! Deallocate old xadj array
    if (allocated(mesh%topo%graph_conn%xadj)) then
      deallocate (mesh%topo%graph_conn%xadj)
    end if

    ! Allocate new adjacency index array xadj based on new vtxdist
    allocate (mesh%topo%graph_conn%xadj(mesh%topo%graph_conn%vtxdist(irank + 2) - mesh%topo%graph_conn%vtxdist(irank + 1) + 1))

    call get_max_faces(max_faces)

    ! Allocate temporary 2D integer work array and initialise to 0
    allocate (tmp_int2d(mesh%topo%graph_conn%vtxdist(irank + 2) - mesh%topo%graph_conn%vtxdist(irank + 1), max_faces + 1))
    tmp_int2d = 0

    ! Allocate array to hold number of neighbours for local cells
    if (allocated(mesh%topo%num_nb)) then
      deallocate (mesh%topo%num_nb)
    end if
    allocate (mesh%topo%num_nb(local_num_cells))


    ! Construct a cell->faces lookup table 
    call get_global_num_cells(global_num_cells)
    call create_shared_array(shared_env, [global_num_cells, max_faces], cell_faces, cell_faces_window)
    if (is_root(shared_env)) then
      allocate (cell_faces_counters(global_num_cells))
      cell_faces(:,:) = -1
      cell_faces_counters(:) = 1
      call get_global_num_faces(global_num_faces)
      do i = 1, global_num_faces
        face_nb1 = mesh%topo%face_cell1(i)
        face_nb2 = mesh%topo%face_cell2(i)

        ! Only add non-boundary cells to the table
        if (face_nb1 /= 0) then
          cell_faces(face_nb1, cell_faces_counters(face_nb1)) = i 
          cell_faces_counters(face_nb1) = cell_faces_counters(face_nb1) + 1
        end if
        if (face_nb2 /= 0) then
          cell_faces(face_nb2, cell_faces_counters(face_nb2)) = i 
          cell_faces_counters(face_nb2) = cell_faces_counters(face_nb2) + 1
        end if
      end do
    end if
    
    call sync(shared_env)
    
    ! Use cell->faces lookup table to compute mesh connectivity for the local cells
    do i = 1, local_num_cells
      global_index_p = mesh%topo%global_indices(i)
      do j = 1, max_faces
        face = cell_faces(global_index_p, j)
        if (face .ne. -1) then

          ! The neighbouring cell is the cell (connected via a face) that is not the same as me
          if (mesh%topo%face_cell1(face) == global_index_p) then
            neighbour = mesh%topo%face_cell2(face)
          else
            neighbour = mesh%topo%face_cell1(face)
          end if

          ! If neighbour is 0, we have a boundary face. Read the boundary id from bnd_rid
          if (neighbour == 0) then
            neighbour = mesh%topo%bnd_rid(face)
          end if
          call compute_connectivity_add_connection(global_index_p, i, neighbour, face, mesh, tmp_int2d)

        end if
      end do
    end do

    ! New number of local connections
    num_connections = sum(tmp_int2d(:, max_faces + 1))
    call dprint("Number of connections after partitioning: " // str(num_connections))

    ! Allocate new adjncy array based on the new number of computed connections
    if (allocated(mesh%topo%graph_conn%adjncy)) then
      deallocate (mesh%topo%graph_conn%adjncy)
    end if

    ! XXX: is adjncy still needed? why is it allocated?
    allocate (mesh%topo%graph_conn%adjncy(num_connections))

    if (allocated(mesh%topo%face_indices)) then
      deallocate (mesh%topo%face_indices)
    end if
    allocate (mesh%topo%face_indices(max_faces, local_num_cells))

    call flatten_connectivity(tmp_int2d, mesh)

    call get_halo_num_cells(halo_num_cells)
    call dprint("Number of halo cells after partitioning: " // str(halo_num_cells))

    call get_total_num_cells(total_num_cells)
    call dprint("Total number of cells (local + halo) after partitioning: " // str(total_num_cells))

    call set_cell_face_indices()

    call set_num_faces(count_mesh_faces())

    call sync(shared_env)
    call destroy_shared_array(shared_env, cell_faces, cell_faces_window)

  end subroutine compute_face_connectivity

  !v Adds a new global index
  !
  !  Reallocates the global index array, updating the total and halo counts in the mesh.
  subroutine add_new_global_index(global_index, mesh)

    integer(ccs_int), intent(in) :: global_index
    type(ccs_mesh), intent(inout) :: mesh

    integer(ccs_int) :: total_num_cells
    integer(ccs_int) :: halo_num_cells

    integer(ccs_int), dimension(:), allocatable :: tmp_global_indices

    call get_total_num_cells(total_num_cells)

    allocate (tmp_global_indices(total_num_cells + 1))
    tmp_global_indices(1:total_num_cells) = mesh%topo%global_indices(1:total_num_cells)
    tmp_global_indices(total_num_cells + 1) = global_index

    ! Update total and halo cell counts
    call set_total_num_cells(total_num_cells + 1)
    call get_total_num_cells(total_num_cells)
    call get_halo_num_cells(halo_num_cells)
    call set_halo_num_cells(halo_num_cells + 1)

    ! Copy extended global indices back into mesh object
    deallocate (mesh%topo%global_indices)
    allocate (mesh%topo%global_indices(total_num_cells))
    mesh%topo%global_indices(:) = tmp_global_indices(:)
    deallocate (tmp_global_indices)

  end subroutine

  module subroutine compute_connectivity_get_local_cells(par_env, mesh)

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the parition

    integer(ccs_int) :: irank
    integer :: i
    integer :: ctr
    integer :: local_num_cells
    integer(ccs_int) :: global_num_cells

    irank = par_env%proc_id

    ! Count the new number of local cells per rank
    local_num_cells = count(mesh%topo%graph_conn%global_partition == irank)

    call set_local_num_cells(local_num_cells)
    call get_local_num_cells(local_num_cells) ! Ensure using value set within mesh
    ! Abort the execution if any rank has 0 local cells
    ! caused by partitioner error
    call dprint("Number of local cells after partitioning: " // str(local_num_cells))
    if (local_num_cells <= 0) then
      call error_abort("ERROR: Zero local cells found.")
    end if

    ! Allocate and then compute global indices
    call get_local_num_cells(local_num_cells)
    if (allocated(mesh%topo%global_indices)) then
      deallocate (mesh%topo%global_indices)
    end if
    allocate (mesh%topo%global_indices(local_num_cells))
    mesh%topo%global_indices(:) = -1 ! This will allow us to check later

    ctr = 1
    associate (irank => par_env%proc_id, &
               partition => mesh%topo%graph_conn%global_partition)
      call get_global_num_cells(global_num_cells)
      do i = 1, global_num_cells
        if (partition(i) == irank) then
          mesh%topo%global_indices(ctr) = i
          ctr = ctr + 1
        end if
      end do
    end associate

    if (ctr /= (local_num_cells + 1)) then
      call error_abort("Didn't find all my cells!")
    end if

    if (minval(mesh%topo%global_indices) < 1) then
      call error_abort("Didn't register all cells properly!")
    end if

    call get_global_num_cells(global_num_cells)
    if (maxval(mesh%topo%global_indices) > global_num_cells) then
      call error_abort("Global index exceeds range!")
    end if

  end subroutine

  subroutine compute_connectivity_add_connection(face_nb1, face_nb1_local_index, face_nb2, face_index, mesh, tmp_int2d)

    integer(ccs_int), intent(in) :: face_nb1             !< Local cell global index
    integer(ccs_int), intent(in) :: face_nb1_local_index !< Local index of face neighbour 1
    integer(ccs_int), intent(in) :: face_nb2             !< Neighbouring cell global index
    integer(ccs_int), intent(in) :: face_index           !< global face index between cell and its neighbour
    type(ccs_mesh), target, intent(inout) :: mesh        !< The mesh for which to compute the partition
    integer, dimension(:, :), intent(inout) :: tmp_int2d !< Temporary connectivity array

    integer :: fctr
    integer(ccs_int) :: max_faces

    call get_max_faces(max_faces)
    fctr = tmp_int2d(face_nb1_local_index, max_faces + 1) + 1 ! Increment number of faces for this cell
    tmp_int2d(face_nb1_local_index, fctr) = face_nb2          ! Store global index of neighbour cell
    tmp_int2d(face_nb1_local_index, max_faces + 1) = fctr     ! Store number of faces for this cell
    mesh%topo%num_nb(face_nb1_local_index) = fctr

    mesh%topo%global_face_indices(fctr, face_nb1) = face_index ! Update face indices to make its order consistent with nb_indices

  end subroutine compute_connectivity_add_connection


  !v Take the 2D connectivity graph and convert to 1D
  !  Note that cell neighbours are still globally numbered at this point.
  subroutine flatten_connectivity(tmp_int2d, mesh)

    use meshing, only: set_halo_num_cells

    integer, dimension(:, :), intent(in) :: tmp_int2d
    type(ccs_mesh), target, intent(inout) :: mesh        !< The mesh for which to compute the partition

    integer :: i, j
    integer, dimension(:), allocatable :: tmp1
    integer, dimension(:), allocatable :: tmp2

    integer :: ctr
    integer, dimension(1) :: local_idx

    integer(ccs_int) :: global_num_cells
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: halo_num_cells
    integer(ccs_int) :: max_faces

    type(cell_locator) :: loc_p
    integer(ccs_int) :: global_index_p

    type(neighbour_locator) :: loc_nb

    call set_halo_num_cells(0)
    call get_halo_num_cells(halo_num_cells)
    call get_local_num_cells(local_num_cells)
    call get_max_faces(max_faces)

    ctr = 1

    allocate (tmp1(local_num_cells))
    if (allocated(mesh%topo%nb_indices)) then
      deallocate (mesh%topo%nb_indices)
    end if
    allocate (mesh%topo%nb_indices(max_faces, local_num_cells))

    tmp1(:) = -1

    ! Initialise neighbour indices
    mesh%topo%nb_indices(:, :) = 0_ccs_int

    call get_halo_num_cells(halo_num_cells)
    call set_total_num_cells(local_num_cells + halo_num_cells)
    do i = 1, local_num_cells
      call create_cell_locator(i, loc_p)

      mesh%topo%graph_conn%xadj(i) = ctr

      ! Loop over connections of cell i
      do j = 1, tmp_int2d(i, max_faces + 1)
        associate (nbidx => tmp_int2d(i, j))
          if ((.not. any(mesh%topo%global_indices == nbidx)) .and. (nbidx .gt. 0)) then
            ! Halo cell
            if (.not. any(tmp1 == nbidx)) then
              ! New halo cell
              ! Copy and extend size of halo cells buffer

              call get_halo_num_cells(halo_num_cells)
              call set_halo_num_cells(halo_num_cells + 1)
              call get_halo_num_cells(halo_num_cells)
              call set_total_num_cells(local_num_cells + halo_num_cells)
              if (halo_num_cells > size(tmp1)) then
                allocate (tmp2(size(tmp1) + local_num_cells))

                tmp2(:) = -1
                tmp2(1:size(tmp1)) = tmp1(1:size(tmp1))

                deallocate (tmp1)
                allocate (tmp1, source=tmp2)
                deallocate (tmp2)
              end if

              tmp1(halo_num_cells) = nbidx
            end if

            local_idx = findloc(tmp1, nbidx)
            call create_neighbour_locator(loc_p, j, loc_nb)
            call set_local_index(local_num_cells + local_idx(1), loc_nb)
          end if

          !local_idx = findloc(tmp1, nbidx)
          !topo%adjncy(ctr) = local_idx(1)

          if (nbidx .lt. 0) then
            ! boundary 'cell'
            call create_neighbour_locator(loc_p, j, loc_nb)
            call set_local_index(nbidx, loc_nb)
          end if

          if (any(mesh%topo%global_indices == nbidx)) then
            ! local in cell
            local_idx = findloc(mesh%topo%global_indices, nbidx)
            call create_neighbour_locator(loc_p, j, loc_nb)
            call set_local_index(local_idx(1), loc_nb)
          end if

          mesh%topo%graph_conn%adjncy(ctr) = nbidx
        end associate

        ctr = ctr + 1
      end do ! End j

    end do
    mesh%topo%graph_conn%xadj(local_num_cells + 1) = ctr

    allocate (tmp2(local_num_cells + halo_num_cells))
    do i = 1, local_num_cells
      call create_cell_locator(i, loc_p)
      call get_global_index(loc_p, global_index_p)
      tmp2(i) = global_index_p
    end do
    do i = 1, halo_num_cells
      tmp2(local_num_cells + i) = tmp1(i)
    end do
    deallocate (mesh%topo%global_indices)
    allocate (mesh%topo%global_indices, source=tmp2)

    deallocate (tmp1)
    deallocate (tmp2)

    call get_global_num_cells(global_num_cells)
    if (minval(mesh%topo%global_indices) < 1) then
      call error_abort("Global index < 0! " // str(minval(mesh%topo%global_indices)))
    end if
    if (maxval(mesh%topo%global_indices) > global_num_cells) then
      call error_abort("Global index > " // str(global_num_cells) // "! " // str(maxval(mesh%topo%global_indices)))
    end if

  end subroutine

  module subroutine compute_partitioner_input_generic(par_env, shared_env, mesh)

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The parallel environment
    type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the parition

    call compute_partitioner_input_generic_topo(par_env, shared_env, mesh%topo)
    
  end subroutine compute_partitioner_input_generic

  subroutine compute_partitioner_input_generic_topo(par_env, shared_env, topo)

    use iso_fortran_env, only: int32

    class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The parallel environment
    type(topology), target, intent(inout) :: topo                              !< The mesh topology for which to compute the parition

    ! Local variables
    integer(ccs_int) :: irank
    integer(ccs_int) :: local_num_cells
    
    call set_topo_object(topo)
    irank = par_env%proc_id

    call compute_partitioner_input_generic_graphconn(par_env, shared_env, topo%global_num_cells, &
                                                     topo%face_cell1, topo%face_cell2, topo%max_faces, &
                                                     topo%graph_conn)

    ! Count and store the number of local cells per rank
    local_num_cells = count(topo%graph_conn%global_partition == irank)
    call set_local_num_cells(local_num_cells)
    call get_local_num_cells(local_num_cells) ! Ensure using true value
    call dprint("Initial number of local cells: " // str(local_num_cells))

    ! Count halo cells
    ! XXX: This will count multiple connections to the same halo cell multiple times, so too would
    !      the original implementation!
    associate(graph_conn => topo%graph_conn)
      topo%halo_num_cells = count(graph_conn%global_partition(graph_conn%adjncy) /= irank)
    end associate
    
    call dprint("Initial number of halo cells: " // str(topo%halo_num_cells))
    topo%total_num_cells = local_num_cells + topo%halo_num_cells
    call dprint("Total number of cells (local + halo): " // str(topo%total_num_cells))

    call nullify_topo_object()
    
  end subroutine compute_partitioner_input_generic_topo

  subroutine compute_partitioner_input_generic_graphconn(par_env, shared_env, global_num_cells_shared_size, &
                                                         face_cell1, face_cell2, max_faces, &
                                                         graph_conn)

    use iso_fortran_env, only: int32

    class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The parallel environment
    integer(ccs_int), intent(in) :: global_num_cells_shared_size               !< The global cell count
    integer(ccs_int), dimension(:), intent(in) :: face_cell1                   !< Face neighbour 1
    integer(ccs_int), dimension(:), intent(in) :: face_cell2                   !< Face neighbour 2
    integer(ccs_int), intent(in) :: max_faces                                  !< Maximum face count per cell
    type(graph_connectivity), target, intent(inout) :: graph_conn              !< The graph connectivity for which to compute the parition

    ! Local variables
    integer(ccs_int), dimension(:, :), allocatable :: tmp_int2d ! Temporary 2D integer array

    integer(ccs_int) :: i, j
    integer(ccs_int) :: irank ! MPI rank ID
    integer(ccs_int) :: isize ! Size of MPI world
    integer(ccs_int) :: start_index
    integer(ccs_int) :: end_index
    integer(ccs_int) :: local_index
    integer(ccs_int) :: num_connections
    integer(ccs_int) :: local_num_cells

    irank = par_env%proc_id
    isize = par_env%num_procs

    start_index = int(graph_conn%vtxdist(irank + 1), int32)
    end_index = int(graph_conn%vtxdist(irank + 2), int32) - 1

    ! Allocate global partition array
    call create_shared_array(shared_env, global_num_cells_shared_size, graph_conn%global_partition, graph_conn%global_partition_window)

    ! Initial global partition
    if (is_root(shared_env)) then
      do i = 1, size(graph_conn%vtxdist) - 1
        j = i - 1
        graph_conn%global_partition(graph_conn%vtxdist(i):graph_conn%vtxdist(i + 1) - 1) = j
      end do
    end if

    call sync(shared_env)

    ! Count the number of local cells per rank
    local_num_cells = count(graph_conn%global_partition == irank)

    ! Allocate adjacency index array xadj based on vtxdist
    allocate (graph_conn%xadj(graph_conn%vtxdist(irank + 2) - graph_conn%vtxdist(irank + 1) + 1))

    ! Allocate temporary 2D integer work array and initialise to 0
    allocate (tmp_int2d(graph_conn%vtxdist(irank + 2) - graph_conn%vtxdist(irank + 1), max_faces + 1))
    tmp_int2d = 0

    call build_connectivity_graph([ start_index, end_index ], face_cell1, face_cell2, tmp_int2d)
    
    num_connections = sum(tmp_int2d(:, max_faces + 1))
    call dprint("Initial number of connections: " // str(num_connections))

    ! Allocate adjncy array based on the number of computed connections
    allocate (graph_conn%adjncy(num_connections))
    ! Allocate local partition array
    allocate (graph_conn%local_partition(graph_conn%vtxdist(irank + 2) - graph_conn%vtxdist(irank + 1)))

    local_index = 1

    do i = 1, end_index - start_index + 1  ! Loop over local cells

      graph_conn%xadj(i) = local_index                          ! Pointer to start of current

      do j = 1, tmp_int2d(i, max_faces + 1)               ! Loop over number of faces
        graph_conn%adjncy(local_index + j - 1) = tmp_int2d(i, j) ! Store global IDs of neighbour cells
      end do

      local_index = local_index + tmp_int2d(i, max_faces + 1)
      graph_conn%xadj(i + 1) = local_index

    end do

    ! Allocate weight arrays
    allocate (graph_conn%adjwgt(num_connections))
    ! Allocate local partition array
    allocate (graph_conn%vwgt(graph_conn%vtxdist(irank + 2) - graph_conn%vtxdist(irank + 1)))

    deallocate (tmp_int2d)

  end subroutine compute_partitioner_input_generic_graphconn

  subroutine build_connectivity_graph(local_range, edge_starts, edge_ends, connectivity_graph)

    integer(ccs_int), dimension(2), intent(in) :: local_range ! Local start/end indices
    integer(ccs_int), dimension(:), intent(in) :: edge_starts ! Starting vertices of each edge
    integer(ccs_int), dimension(:), intent(in) :: edge_ends   ! Ending vertices of each edge
    integer(ccs_int), dimension(:, :), intent(out) :: connectivity_graph

    ! Local variables
    integer(ccs_int) :: max_degree ! Maximum number of edges connected to a vertex
    integer(ccs_int) :: nedges ! Number of edges in the graph
    integer(ccs_int) :: start_index ! Start index of local range
    integer(ccs_int) :: end_index   ! End index of local range

    integer(ccs_int) :: vtx1, vtx2
    
    integer(ccs_int) :: local_index
    integer(ccs_int) :: i
    integer(ccs_int) :: k
    
    max_degree = size(connectivity_graph, 2) - 1

    nedges = size(edge_starts)
    if (size(edge_ends) /= nedges) then
       call error_abort("Size of edge start/end arrays doesn't match")
    end if

    start_index = local_range(1)
    end_index = local_range(2)
    
    ! All ranks loop over all the graph edges
    do i = 1, nedges

      vtx1 = edge_starts(i)
      vtx2 = edge_ends(i)

      if (vtx2 /= 0) then
        if (vtx1 /= 0) then
          ! If edge vertex 1 is local to the current rank
          if (vtx1 >= start_index) then
            if (vtx1 <= end_index) then
              local_index = vtx1 - start_index + 1                    ! Local vertex index
              k = connectivity_graph(local_index, max_degree + 1) + 1 ! Increment number of edges
                                                                      ! for this vertex
              connectivity_graph(local_index, k) = vtx2               ! Store global index of neighbour
                                                                      ! vertex
              connectivity_graph(local_index, max_degree + 1) = k     ! Store number of edges for this
                                                                      ! vertex
            end if
          end if

          ! If edge vertex 2 is local to the current rank
          if (vtx2 >= start_index) then
            if (vtx2 <= end_index) then
              local_index = vtx2 - start_index + 1                    ! Local vertex index
              k = connectivity_graph(local_index, max_degree + 1) + 1 ! Increment number of edges
                                                                      ! for this vertex
              connectivity_graph(local_index, k) = vtx1               ! Store global index of neighbour
                                                                      ! vertex
              connectivity_graph(local_index, max_degree + 1) = k     ! Store number of edges for this
                                                                      ! vertex
            end if
          end if
        end if
      else
        ! If edge vertex 2 is 0 we have a boundary edge
      end if

    end do
    
  end subroutine build_connectivity_graph
  
  !v Deallocate partitioner data structures associated with the mesh
  module subroutine cleanup_partitioner_data(shared_env, mesh)
    
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The parallel environment
    type(ccs_mesh), target, intent(inout) :: mesh   !< The mesh

    call cleanup_partitioner_data_topo(shared_env, mesh%topo)

  end subroutine cleanup_partitioner_data
  subroutine cleanup_partitioner_data_topo(shared_env, topo)
    
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The parallel environment
    type(topology), target, intent(inout) :: topo   !< The mesh topology

    call cleanup_partitioner_data_graphconn(shared_env, topo%graph_conn)

  end subroutine cleanup_partitioner_data_topo
  subroutine cleanup_partitioner_data_graphconn(shared_env, graph_conn)
    
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The parallel environment
    type(graph_connectivity), target, intent(inout) :: graph_conn   !< The mesh topology graph connectivity

    if (allocated(graph_conn%vtxdist)) then 
      deallocate(graph_conn%vtxdist)
      call dprint("graph_conn%vtxdist deallocated.")
    end if

    if (allocated(graph_conn%xadj)) then 
      deallocate(graph_conn%xadj)
      call dprint("graph_conn%xadj deallocated.")
    end if
    
    if (allocated(graph_conn%adjncy)) then 
      deallocate(graph_conn%adjncy)
      call dprint("graph_conn%adjncy deallocated.")
    end if

    if (allocated(graph_conn%local_partition)) then 
      deallocate(graph_conn%local_partition)
      call dprint("graph_conn%local_partition deallocated.")
    end if

    if (allocated(graph_conn%adjwgt)) then 
      deallocate(graph_conn%adjwgt)
      call dprint("graph_conn%adjwgt deallocated.")
    end if

    if (allocated(graph_conn%vwgt)) then 
      deallocate(graph_conn%vwgt)
      call dprint("topo%vwgt deallocated.")
    end if

    if (associated(graph_conn%global_partition)) then 
      call destroy_shared_array(shared_env, graph_conn%global_partition, graph_conn%global_partition_window)
      call dprint("graph_conn%global_partition deallocated.")
    end if

  end subroutine cleanup_partitioner_data_graphconn

  !v Compute and report the partitioning quality.
  !
  !  The following metrics are implemented
  !  - The "surface to volume ratio" nhalo / nlocal (averaged)
  !  - The minimum departure from load balance min(nlocal) / avg(nlocal)
  !  - The maximum departure from load balance max(nlocal) / avg(nlocal)
  module subroutine print_partition_quality(par_env)

    use case_config, only : compute_partqual
    
    class(parallel_environment), intent(in) :: par_env

    real(ccs_real) :: s2v ! Surface to volume ratio
    real(ccs_real) :: ulb ! Under load balance (minimum)
    real(ccs_real) :: olb ! Over load balance (maximum)

    if (compute_partqual) then
      call compute_partition_quality(par_env, s2v, ulb, olb)
      if (is_root(par_env)) then
        print *, "Partitioning report:"
        print *, "- Surface:Volume ratio:", s2v
        print *, "- Under load balance:", ulb
        print *, "- Over load balance:", olb
      end if
    end if
    
  end subroutine print_partition_quality

  !v Compute the partitioning quality.
  !
  !  The following metrics are implemented
  !  - The "surface to volume ratio" nhalo / nlocal (averaged)
  !  - The minimum departure from load balance min(nlocal) / avg(nlocal)
  !  - The maximum departure from load balance max(nlocal) / avg(nlocal)
  subroutine compute_partition_quality(par_env, s2v, ulb, olb)

    use mpi
    
    class(parallel_environment), intent(in) :: par_env
    real(ccs_real), intent(out) :: s2v ! Surface to volume ratio
    real(ccs_real), intent(out) :: ulb ! Under load balance (minimum)
    real(ccs_real), intent(out) :: olb ! Over load balance (maximum)

    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: halo_num_cells

    integer(ccs_int) :: local_num_cells_stat
    
    real(ccs_real) :: local_num_cells_avg

    integer(ccs_err) :: ierr

    call get_local_num_cells(local_num_cells)
    call get_halo_num_cells(halo_num_cells)

    associate(nprocs => real(par_env%num_procs, ccs_real))
      ! Compute average surface to volume ratio
      s2v = real(halo_num_cells, ccs_real) / real(local_num_cells, ccs_real)
      select type(par_env)
      type is (parallel_environment_mpi)
        call MPI_Allreduce(MPI_IN_PLACE, s2v, 1, MPI_DOUBLE_PRECISION, MPI_SUM, par_env%comm, ierr)
      class default
        call error_abort("Unsupported parallel environment")
      end select
      s2v = s2v / nprocs

      ! Compute average local cell count
      select type(par_env)
      type is (parallel_environment_mpi)
        call MPI_Allreduce(local_num_cells, local_num_cells_stat, 1, MPI_INTEGER, MPI_SUM, par_env%comm, ierr)
      class default
        call error_abort("Unsupported parallel environment")
      end select
      local_num_cells_avg = real(local_num_cells_stat, ccs_real) / nprocs

      ! Compute under load balance
      select type(par_env)
      type is (parallel_environment_mpi)
        call MPI_Allreduce(local_num_cells, local_num_cells_stat, 1, MPI_INTEGER, MPI_MIN, par_env%comm, ierr)
      class default
        call error_abort("Unsupported parallel environment")
      end select
      ulb = real(local_num_cells_stat, ccs_real) / local_num_cells_avg

      ! Compute over load balance
      select type(par_env)
      type is (parallel_environment_mpi)
        call MPI_Allreduce(local_num_cells, local_num_cells_stat, 1, MPI_INTEGER, MPI_MAX, par_env%comm, ierr)
      class default
        call error_abort("Unsupported parallel environment")
      end select
      olb = real(local_num_cells_stat, ccs_real) / local_num_cells_avg
    end associate
    
  end subroutine compute_partition_quality
  
end submodule
