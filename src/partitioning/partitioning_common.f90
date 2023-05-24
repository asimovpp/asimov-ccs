submodule(partitioning) partitioning_common
#include "ccs_macros.inc"

  use kinds, only: ccs_int, ccs_err
  use types, only: cell_locator, neighbour_locator, vertex_neighbour_locator
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
                     get_count_vertex_neighbours

  implicit none

contains

  ! Compute the new topology connectivity after partitioning
  module subroutine compute_connectivity(par_env, mesh)

    use mpi

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the parition

    ! Local variables
    integer(ccs_int) :: irank ! MPI rank ID
    integer(ccs_int) :: isize ! Size of MPI world
    integer(ccs_int) :: i
    
    irank = par_env%proc_id
    isize = par_env%num_procs

    ! Prepare global data needed later
    call store_global_vertex_connectivity(par_env, mesh)
    

    ! Get global indices of local cells
    call compute_connectivity_get_local_cells(par_env, mesh)

    ! Recompute vtxdist array based on the new partition
    mesh%topo%vtxdist(1) = 1
    do i = 2, isize + 1
      mesh%topo%vtxdist(i) = count(mesh%topo%global_partition == (i - 2)) + mesh%topo%vtxdist(i - 1)
    end do

    if (irank == 0) then
      do i = 1, isize + 1
        call dprint("new vtxdist(" // str(i) // "): " // str(int(mesh%topo%vtxdist(i))))
      end do
    end if

    call compute_face_connectivity(par_env, mesh)
    call compute_vertex_connectivity(par_env, mesh)
    
  end subroutine compute_connectivity

  subroutine compute_face_connectivity(par_env, mesh)

    use iso_fortran_env, only: int32

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(ccs_mesh), target, intent(inout) :: mesh                           !< The mesh for which to compute the parition

    ! Local variables
    integer(ccs_int), dimension(:, :), allocatable :: tmp_int2d ! Temporary 2D integer array
    integer(ccs_int) :: irank ! MPI rank ID
    integer(ccs_int) :: isize ! Size of MPI world
    integer(ccs_int) :: i
    integer(ccs_int) :: start_index
    integer(ccs_int) :: end_index
    integer(ccs_int) :: face_nb1
    integer(ccs_int) :: face_nb2
    integer(ccs_int) :: num_connections
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: global_num_cells
    integer(ccs_int) :: halo_num_cells
    integer(ccs_int) :: total_num_cells
    integer(ccs_int) :: global_num_faces
    integer(ccs_int) :: max_faces
    
    irank = par_env%proc_id
    isize = par_env%num_procs

    call get_local_num_cells(mesh, local_num_cells)
    
    ! Deallocate old xadj array
    if (allocated(mesh%topo%xadj)) then
      deallocate (mesh%topo%xadj)
    end if

    ! Allocate new adjacency index array xadj based on new vtxdist
    allocate (mesh%topo%xadj(mesh%topo%vtxdist(irank + 2) - mesh%topo%vtxdist(irank + 1) + 1))

    if (allocated(mesh%topo%global_boundaries) .eqv. .false.) then
      call get_global_num_cells(mesh, global_num_cells)
      allocate (mesh%topo%global_boundaries(global_num_cells))
    end if

    ! Reset global_boundaries array
    mesh%topo%global_boundaries = 0

    call get_max_faces(mesh, max_faces)

    ! Allocate temporary 2D integer work array and initialise to 0
    allocate (tmp_int2d(mesh%topo%vtxdist(irank + 2) - mesh%topo%vtxdist(irank + 1), max_faces + 1))
    tmp_int2d = 0

    start_index = int(mesh%topo%vtxdist(irank + 1), int32)
    end_index = int(mesh%topo%vtxdist(irank + 2), int32) - 1

    ! Allocate array to hold number of neighbours for local cells
    if (allocated(mesh%topo%num_nb)) then
      deallocate (mesh%topo%num_nb)
    end if
    allocate (mesh%topo%num_nb(local_num_cells))

    ! All ranks loop over all the faces again
    call get_global_num_faces(mesh, global_num_faces)
    do i = 1, global_num_faces

      face_nb1 = mesh%topo%face_cell1(i)
      face_nb2 = mesh%topo%face_cell2(i)

      ! If face neighbour 1 is local to the current rank
      ! and face neighbour 2 is not 0
      if (any(mesh%topo%global_indices == face_nb1) .and. (face_nb2 .ne. 0)) then
        call compute_connectivity_add_connection(face_nb1, face_nb2, mesh, tmp_int2d)
      end if

      ! If face neighbour 2 is local to the current rank
      ! and face neighbour 1 is not 0
      if (any(mesh%topo%global_indices == face_nb2) .and. (face_nb1 .ne. 0)) then
        call compute_connectivity_add_connection(face_nb2, face_nb1, mesh, tmp_int2d)
      end if

      ! If face neighbour 1 is local and if face neighbour 2 is 0 we have a boundary face
      if (any(mesh%topo%global_indices == face_nb1) .and. (face_nb2 .eq. 0)) then
        mesh%topo%global_boundaries(face_nb1) = mesh%topo%global_boundaries(face_nb1) + 1

        ! read the boundary id from bnd_rid
        face_nb2 = mesh%topo%bnd_rid(i)
        call compute_connectivity_add_connection(face_nb1, face_nb2, mesh, tmp_int2d)
      end if

    end do

    ! New number of local connections
    num_connections = sum(tmp_int2d(:, max_faces + 1))
    call dprint("Number of connections after partitioning: " // str(num_connections))

    ! Allocate new adjncy array based on the new number of computed connections
    if (allocated(mesh%topo%adjncy)) then
      deallocate (mesh%topo%adjncy)
    end if

    allocate (mesh%topo%adjncy(num_connections))

    if (allocated(mesh%topo%face_indices)) then
      deallocate (mesh%topo%face_indices)
    end if
    allocate (mesh%topo%face_indices(max_faces, local_num_cells))

    call flatten_connectivity(tmp_int2d, mesh)

    call get_halo_num_cells(mesh, halo_num_cells)
    call dprint("Number of halo cells after partitioning: " // str(halo_num_cells))

    call get_total_num_cells(mesh, total_num_cells)
    call dprint("Total number of cells (local + halo) after partitioning: " // str(total_num_cells))

    call set_cell_face_indices(mesh)

    call set_num_faces(count_mesh_faces(mesh), mesh)
    
  end subroutine compute_face_connectivity

  subroutine store_global_vertex_connectivity(par_env, mesh)

    use mpi
    
    class(parallel_environment), intent(in) :: par_env
    type(ccs_mesh), intent(inout) :: mesh

    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: global_num_cells
    integer(ccs_int) :: max_vert_nb

    integer(ccs_int), dimension(:), allocatable :: tmp_arr
    integer(ccs_int) :: i, j
    integer(ccs_int) :: idx
    integer(ccs_int) :: vert_nb_idx
    
    integer(ccs_err) :: ierr
    
    global_num_cells = mesh%topo%global_num_cells

    if (size(mesh%topo%vert_nb_indices, 2) /= global_num_cells) then
      ! We only have a local view of vertex neighbours

      max_vert_nb = maxval(mesh%topo%num_vert_nb)
      select type(par_env)
      type is(parallel_environment_mpi)
        call MPI_Allreduce(MPI_IN_PLACE, max_vert_nb, 1, MPI_INTEGER, MPI_MAX, par_env%comm, ierr)
      class default
        call error_abort("Unsupported parallel environment!")
      end select

      allocate(tmp_arr(max_vert_nb * global_num_cells))
      tmp_arr(:) = 0

      ! XXX: cannot read local_num_cells - arrays haven't been resized!
      local_num_cells = size(mesh%topo%num_vert_nb)

      do i = 1, local_num_cells
        idx = max_vert_nb * (mesh%topo%global_indices(i) - 1)
        
        do j = 1, mesh%topo%num_vert_nb(i)
          vert_nb_idx = mesh%topo%vert_nb_indices(j, i)
          if (vert_nb_idx > 0) then
            tmp_arr(idx + j) = mesh%topo%global_indices(vert_nb_idx)
          else
            tmp_arr(idx + j) = vert_nb_idx
          end if
        end do
      end do

      select type(par_env)
      type is(parallel_environment_mpi)
        call MPI_Allreduce(MPI_IN_PLACE, tmp_arr, max_vert_nb * global_num_cells, MPI_INTEGER, MPI_SUM, par_env%comm, ierr)
      class default
        call error_abort("Unsupported parallel environment!")
      end select

      if (allocated(mesh%topo%vert_nb_indices)) then
        deallocate(mesh%topo%vert_nb_indices)
      end if
      allocate(mesh%topo%vert_nb_indices(max_vert_nb, global_num_cells))
      mesh%topo%vert_nb_indices(:, :) = 0
      do i = 1, global_num_cells
        do j = 1, max_vert_nb
          idx = max_vert_nb * (i - 1) + j
          mesh%topo%vert_nb_indices(j, i) = tmp_arr(idx)
        end do
      end do
      
      deallocate(tmp_arr)
   else
      if (par_env%proc_id == par_env%root) then
         print *, "Continuing without constructing global vertex neighbours"
      end if
   end if
    
  end subroutine store_global_vertex_connectivity
  
  subroutine compute_vertex_connectivity(par_env, mesh)

    use mpi

    class(parallel_environment), intent(in) :: par_env
    type(ccs_mesh), target, intent(inout) :: mesh !< The mesh for which to compute the parition

    integer(ccs_int) :: local_num_cells

    integer(ccs_int) :: i, j
    integer(ccs_int) :: local_idx
    integer(ccs_int) :: global_idx
    
    integer(ccs_int) :: max_vert_nb
    
    integer(ccs_err) :: ierr

    integer(ccs_int), dimension(:, :), allocatable :: tmp_2d
    integer(ccs_int), dimension(:), allocatable :: tmp_1d

    type(cell_locator) :: loc_p
    integer(ccs_int) :: nvnb

    integer(ccs_int) :: total_num_cells
    integer(ccs_int) :: halo_num_cells

    integer(ccs_int) :: global_num_cells
    
    if (.not. allocated(mesh%topo%global_vertex_indices)) then
      call error_abort("The global vertex indices array was deallocated prematurely!")
    else
      if (size(mesh%topo%global_vertex_indices, 2) /= mesh%topo%global_num_cells) then
        call error_abort("Need the global cell-vertex connnectivity, not just local!")
      end if
    end if

    call get_global_num_cells(mesh, global_num_cells)
    call get_local_num_cells(mesh, local_num_cells)

    max_vert_nb = maxval(mesh%topo%num_vert_nb)
    select type(par_env)
    type is(parallel_environment_mpi)
      call MPI_Allreduce(MPI_IN_PLACE, max_vert_nb, 1, MPI_INTEGER, MPI_MAX, par_env%comm, ierr)
    class default
      call error_abort("Unsupported parallel environment!")
    end select

    !! XXX: Need to get the maximum number of vertex neighbours BEFORE deallocating
    if (allocated(mesh%topo%num_vert_nb)) then
      deallocate(mesh%topo%num_vert_nb)
    end if
    allocate(mesh%topo%num_vert_nb(local_num_cells))

    ! Check vertex neighbours as input
    if (any(mesh%topo%vert_nb_indices > global_num_cells)) then
       call error_abort("Global vertex neighbour indices > global_num_cells")
    end if
    
    ! Copy vertex neighbour indices from global array
    allocate(tmp_2d, source=mesh%topo%vert_nb_indices)
    deallocate(mesh%topo%vert_nb_indices)
    allocate(mesh%topo%vert_nb_indices(max_vert_nb, local_num_cells))
    do local_idx = 1, local_num_cells
      call create_cell_locator(mesh, local_idx, loc_p)
      call get_global_index(loc_p, global_idx)

      mesh%topo%vert_nb_indices(:, local_idx) = tmp_2d(:, global_idx)
      mesh%topo%num_vert_nb(local_idx) = count(mesh%topo%vert_nb_indices(:, local_idx) /= 0)
    end do
    deallocate(tmp_2d)
    
    ! Convert global->local indices
    do i = 1, local_num_cells
      call create_cell_locator(mesh, i, loc_p)
      
      !call get_count_vertex_neighbours(loc_p, nvnb)
      nvnb = mesh%topo%num_vert_nb(i)
      do j = 1, nvnb
        ! We can't use the neighbour index because currently we have global indices and this breaks the error checking.
        global_idx = mesh%topo%vert_nb_indices(j, i)
        
        if (global_idx > 0) then
          local_idx = findloc(mesh%topo%global_indices, global_idx, 1)
          if (local_idx == 0) then
            ! New global index
            call get_total_num_cells(mesh, total_num_cells)

            allocate(tmp_1d(total_num_cells + 1))
            tmp_1d(1:total_num_cells) = mesh%topo%global_indices(1:total_num_cells)
            tmp_1d(total_num_cells + 1) = global_idx

            ! Update total and halo cell counts
            call set_total_num_cells(total_num_cells + 1, mesh)
            call get_total_num_cells(mesh, total_num_cells)
            call get_halo_num_cells(mesh, halo_num_cells)
            call set_halo_num_cells(halo_num_cells + 1, mesh)

            ! Copy extended global indices back into mesh object
            deallocate(mesh%topo%global_indices)
            allocate(mesh%topo%global_indices(total_num_cells))
            mesh%topo%global_indices(:) = tmp_1d(:)
            deallocate(tmp_1d)

            local_idx = total_num_cells
          end if

          ! As above, can't use vertex neighbour locators as currently have global indices in the arrays.
          mesh%topo%vert_nb_indices(j, i) = local_idx
        end if
      end do
    end do
    
  end subroutine compute_vertex_connectivity
  
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
    local_num_cells = count(mesh%topo%global_partition == irank)
    call set_local_num_cells(local_num_cells, mesh)
    call get_local_num_cells(mesh, local_num_cells) ! Ensure using value set within mesh
    call dprint("Number of local cells after partitioning: " // str(local_num_cells))

    ! Allocate and then compute global indices
    if (allocated(mesh%topo%global_indices)) then
      deallocate (mesh%topo%global_indices)
    end if
    call get_local_num_cells(mesh, local_num_cells)
    if (allocated(mesh%topo%global_indices)) then
       deallocate(mesh%topo%global_indices)
    end if
    allocate (mesh%topo%global_indices(local_num_cells))
    mesh%topo%global_indices(:) = -1 ! This will allow us to check later

    ctr = 1
    associate (irank => par_env%proc_id, &
               partition => mesh%topo%global_partition)
      call get_global_num_cells(mesh, global_num_cells)
      do i = 1, global_num_cells
        if (partition(i) == irank) then
          mesh%topo%global_indices(ctr) = i
          ctr = ctr + 1
        end if
      end do
    end associate

    if (ctr /= (local_num_cells + 1)) then
      print *, "ERROR: didn't find all my cells!"
      stop
    end if

    if (minval(mesh%topo%global_indices) < 1) then
      print *, "ERROR: didn't register all cells properly!"
      stop
    end if

    call get_global_num_cells(mesh, global_num_cells)
    if (maxval(mesh%topo%global_indices) > global_num_cells) then
      print *, "ERROR: global index exceeds range!"
      stop
    end if

  end subroutine

  subroutine compute_connectivity_add_connection(face_nb1, face_nb2, mesh, tmp_int2d)

    integer(ccs_int), intent(in) :: face_nb1             !< Local cell global index
    integer(ccs_int), intent(in) :: face_nb2             !< Neighbouring cell global index
    type(ccs_mesh), target, intent(inout) :: mesh        !< The mesh for which to compute the partition
    integer, dimension(:, :), intent(inout) :: tmp_int2d !< Temporary connectivity array

    integer, dimension(1) :: local_index
    integer :: fctr

    integer(ccs_int) :: max_faces

    local_index = findloc(mesh%topo%global_indices, face_nb1)
    if (local_index(1) <= 0) then
      print *, "ERROR: failed to find face neighbour in global indices, findloc: "
      print *, "- ANY: ", any(mesh%topo%global_indices == face_nb1)
      print *, "- local_index: ", local_index
      print *, "- Face neighbour index: ", face_nb1
      print *, "- Global indices: ", mesh%topo%global_indices
      stop
    end if

    call get_max_faces(mesh, max_faces)
    fctr = tmp_int2d(local_index(1), max_faces + 1) + 1 ! Increment number of faces for this cell
    tmp_int2d(local_index(1), fctr) = face_nb2               ! Store global index of neighbour cell
    tmp_int2d(local_index(1), max_faces + 1) = fctr     ! Store number of faces for this cell
    mesh%topo%num_nb(local_index(1)) = fctr

  end subroutine

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

    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: halo_num_cells
    integer(ccs_int) :: max_faces

    type(cell_locator) :: loc_p
    integer(ccs_int) :: global_index_p

    type(neighbour_locator) :: loc_nb

    call set_halo_num_cells(0, mesh)
    call get_halo_num_cells(mesh, halo_num_cells)
    call get_local_num_cells(mesh, local_num_cells)
    call get_max_faces(mesh, max_faces)

    ctr = 1

    allocate (tmp1(local_num_cells))
    if (allocated(mesh%topo%nb_indices)) then
      deallocate (mesh%topo%nb_indices)
    end if
    allocate (mesh%topo%nb_indices(max_faces, local_num_cells))

    tmp1(:) = -1

    ! Initialise neighbour indices
    mesh%topo%nb_indices(:, :) = 0_ccs_int

    call get_halo_num_cells(mesh, halo_num_cells)
    call set_total_num_cells(local_num_cells + halo_num_cells, mesh)
    do i = 1, local_num_cells
      call create_cell_locator(mesh, i, loc_p)

      mesh%topo%xadj(i) = ctr

      ! Loop over connections of cell i
      do j = 1, tmp_int2d(i, max_faces + 1)
        associate (nbidx => tmp_int2d(i, j))
          if ((.not. any(mesh%topo%global_indices == nbidx)) .and. (nbidx .gt. 0)) then
            ! Halo cell
            if (.not. any(tmp1 == nbidx)) then
              ! New halo cell
              ! Copy and extend size of halo cells buffer

              call get_halo_num_cells(mesh, halo_num_cells)
              call set_halo_num_cells(halo_num_cells + 1, mesh)
              call get_halo_num_cells(mesh, halo_num_cells)
              call set_total_num_cells(local_num_cells + halo_num_cells, mesh)
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

          mesh%topo%adjncy(ctr) = nbidx
        end associate

        ctr = ctr + 1
      end do ! End j

    end do
    mesh%topo%xadj(local_num_cells + 1) = ctr

    allocate (tmp2(local_num_cells + halo_num_cells))
    do i = 1, local_num_cells
      call create_cell_locator(mesh, i, loc_p)
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

  end subroutine

end submodule
