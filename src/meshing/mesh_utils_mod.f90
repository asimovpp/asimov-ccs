module mesh_utils
#include "ccs_macros.inc"

  use constants, only: ndim

  use utils, only: exit_print
  use kinds, only: ccs_int, ccs_real
  use types, only: ccs_mesh, face_locator, cell_locator, neighbour_locator
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi
  use meshing, only: get_global_index, get_local_index, count_neighbours, &
                     set_cell_location, set_neighbour_location, set_face_location, &
                     set_face_index, get_boundary_status, get_local_status
  use bc_constants

  implicit none

  !^ @note Named constants for faces of hexahedral cells follow the convention that the lower
  !        boundary on a given axis is numbered first, i.e.
  !                                
  !         +----------+
  !        /|    4    /| 
  !       +----------+ |
  !       | |        | |
  !     1 | |        | | 2
  !       | +--------|-+
  !       |/    3    |/
  !       +----------+
  !             
  !
  !  @endnote
  integer, parameter :: left = 1_ccs_int
  integer, parameter :: right = 2_ccs_int
  integer, parameter :: bottom = 3_ccs_int
  integer, parameter :: top = 4_ccs_int
  integer, parameter :: back = 5_ccs_int
  integer, parameter :: front = 6_ccs_int


  private
  public :: build_square_mesh
  public :: build_mesh
  public :: global_start
  public :: local_count
  public :: count_mesh_faces

contains

  !v Utility constructor to build a square mesh.
  !
  !  Builds a Cartesian grid of NxN cells on the domain LxL.
  function build_square_mesh(par_env, cps, side_length) result(mesh)

    class(parallel_environment), intent(in) :: par_env !< The parallel environment to construct the mesh.
    integer(ccs_int), intent(in) :: cps                !< Number of cells per side of the mesh.
    real(ccs_real), intent(in) :: side_length          !< The length of each side.

    type(ccs_mesh) :: mesh                             !< The resulting mesh.

    integer(ccs_int) :: start_global    ! The (global) starting index of a partition
    integer(ccs_int) :: end_global      ! The (global) last index of a partition
    integer(ccs_int) :: i               ! Loop counter
    integer(ccs_int) :: ii              ! Zero-indexed loop counter (simplifies some operations)
    integer(ccs_int) :: index_counter   ! Local index counter
    integer(ccs_int) :: face_counter    ! Cell-local face counter
    integer(ccs_int) :: comm_rank       ! The process ID within the parallel environment
    integer(ccs_int) :: comm_size       ! The size of the parallel environment

    integer(ccs_int) :: index_nb        ! The local index of a neighbour cell
    integer(ccs_int) :: global_index_nb ! The global index of a neighbour cell

    select type (par_env)
    type is (parallel_environment_mpi)

      ! Set the global mesh parameters
      mesh%topo%global_num_cells = cps**2
      mesh%geo%h = side_length / real(cps, ccs_real)

      ! Associate aliases to make code easier to read
      associate (nglobal => mesh%topo%global_num_cells, &
                 h => mesh%geo%h)

        ! Determine ownership range
        comm_rank = par_env%proc_id
        comm_size = par_env%num_procs
        start_global = global_start(nglobal, par_env%proc_id, par_env%num_procs)
        mesh%topo%local_num_cells = local_count(nglobal, par_env%proc_id, par_env%num_procs)
        end_global = start_global + (mesh%topo%local_num_cells - 1)

        ! Set max faces per cell (constant, 4)
        mesh%topo%max_faces = 4_ccs_int
        
        ! Allocate mesh arrays
        allocate (mesh%topo%global_indices(mesh%topo%local_num_cells))
        allocate (mesh%topo%num_nb(mesh%topo%local_num_cells))
        allocate (mesh%topo%nb_indices(mesh%topo%max_faces, mesh%topo%local_num_cells))
        allocate (mesh%topo%face_indices(mesh%topo%max_faces, mesh%topo%local_num_cells))

        ! Initialise mesh arrays
    
        mesh%topo%num_nb(:) = mesh%topo%max_faces ! All cells have 4 neighbours (possibly ghost/boundary cells)

        ! First set the global index of local cells
        index_counter = 1_ccs_int
        do i = start_global, end_global
          mesh%topo%global_indices(index_counter) = i
          index_counter = index_counter + 1
        end do

        ! Assemble cells and faces
        ! XXX: Negative neighbour indices are used to indicate boundaries using the same numbering
        !      as cell-relative neighbour indexing, i.e.
        !        -1 = left boundary
        !        -2 = right boundary
        !        -3 = bottom boundary
        !        -4 = top boundary
        index_counter = 1_ccs_int ! Set local indexing starting from 1...n
        do i = start_global, end_global
          ii = i - 1_ccs_int

          ! Construct left (1) face/neighbour
          face_counter = left
          if (modulo(ii, cps) == 0_ccs_int) then
            index_nb = -left
            global_index_nb = -left
          else
            index_nb = index_counter - 1_ccs_int
            global_index_nb = i - 1_ccs_int
          end if
          call build_local_mesh_add_neighbour(index_counter, face_counter, index_nb, global_index_nb, mesh)

          ! Construct right (2) face/neighbour
          face_counter = right
          if (modulo(ii, cps) == (cps - 1_ccs_int)) then
            index_nb = -right
            global_index_nb = -right
          else
            index_nb = index_counter + 1_ccs_int
            global_index_nb = i + 1_ccs_int
          end if
          call build_local_mesh_add_neighbour(index_counter, face_counter, index_nb, global_index_nb, mesh)

          ! Construct bottom (3) face/neighbour
          face_counter = bottom
          if ((ii / cps) == 0_ccs_int) then
            index_nb = -bottom
            global_index_nb = -bottom
          else
            index_nb = index_counter - cps
            global_index_nb = i - cps
          end if
          call build_local_mesh_add_neighbour(index_counter, face_counter, index_nb, global_index_nb, mesh)

          ! Construct top (4) face/neighbour
          face_counter = top
          if ((ii / cps) == (cps - 1_ccs_int)) then
            index_nb = -top
            global_index_nb = -top
          else
            index_nb = index_counter + cps
            global_index_nb = i + cps
          end if
          call build_local_mesh_add_neighbour(index_counter, face_counter, index_nb, global_index_nb, mesh)

          index_counter = index_counter + 1_ccs_int
        end do
      end associate

      mesh%topo%total_num_cells = size(mesh%topo%global_indices)
      mesh%topo%halo_num_cells = mesh%topo%total_num_cells - mesh%topo%local_num_cells

      allocate (mesh%geo%x_p(ndim, mesh%topo%total_num_cells))
      allocate (mesh%geo%x_f(ndim, mesh%topo%max_faces, mesh%topo%local_num_cells)) !< @note Currently hardcoded as a 2D mesh. @endnote
      allocate (mesh%geo%volumes(mesh%topo%total_num_cells))
      allocate (mesh%geo%face_areas(mesh%topo%max_faces, mesh%topo%local_num_cells))
      allocate (mesh%geo%face_normals(ndim, mesh%topo%max_faces, mesh%topo%local_num_cells)) ! Currently hardcoded as a 2D mesh.

      mesh%geo%volumes(:) = mesh%geo%h**2 !< @note Mesh is square and 2D @endnote
      mesh%geo%face_normals(:, :, :) = 0.0_ccs_real
      mesh%geo%x_p(:, :) = 0.0_ccs_real
      mesh%geo%x_f(:, :, :) = 0.0_ccs_real
      mesh%geo%face_areas(:, :) = mesh%geo%h  ! Mesh is square and 2D

      associate (h => mesh%geo%h)
        do i = 1_ccs_int, mesh%topo%total_num_cells
          ii = mesh%topo%global_indices(i)

          associate (x_p => mesh%geo%x_p(:, i))
            ! Set cell centre
            x_p(1) = (modulo(ii - 1, cps) + 0.5_ccs_real) * h
            x_p(2) = ((ii - 1) / cps + 0.5_ccs_real) * h
          end associate
        end do

        do i = 1_ccs_int, mesh%topo%local_num_cells
          associate (x_p => mesh%geo%x_p(:, i), &
                     x_f => mesh%geo%x_f(:, :, i), &
                     normal => mesh%geo%face_normals(:, :, i))

            face_counter = left
            x_f(1, face_counter) = x_p(1) - 0.5_ccs_real * h
            x_f(2, face_counter) = x_p(2)
            normal(1, face_counter) = -1.0_ccs_real
            normal(2, face_counter) = 0.0_ccs_real

            face_counter = right
            x_f(1, face_counter) = x_p(1) + 0.5_ccs_real * h
            x_f(2, face_counter) = x_p(2)
            normal(1, face_counter) = 1.0_ccs_real
            normal(2, face_counter) = 0.0_ccs_real

            face_counter = bottom
            x_f(1, face_counter) = x_p(1)
            x_f(2, face_counter) = x_p(2) - 0.5_ccs_real * h
            normal(1, face_counter) = 0.0_ccs_real
            normal(2, face_counter) = -1.0_ccs_real

            face_counter = top
            x_f(1, face_counter) = x_p(1)
            x_f(2, face_counter) = x_p(2) + 0.5_ccs_real * h
            normal(1, face_counter) = 0.0_ccs_real
            normal(2, face_counter) = 1.0_ccs_real
          end associate
        end do
      end associate

      mesh%topo%num_faces = count_mesh_faces(mesh)

      call set_cell_face_indices(mesh)

    class default
      call error_abort("Unknown parallel environment type.")

    end select
  end function build_square_mesh


  !v Utility constructor to build a 3D mesh with hex cells.
  !
  !  Builds a Cartesian grid of nx*ny*nz cells.
  function build_mesh(par_env, nx, ny, nz, side_length) result(mesh)

    class(parallel_environment), intent(in) :: par_env !< The parallel environment to construct the mesh.
    integer(ccs_int), intent(in) :: nx                 !< Number of cells in the x direction.
    integer(ccs_int), intent(in) :: ny                 !< Number of cells in the y direction.
    integer(ccs_int), intent(in) :: nz                 !< Number of cells in the z direction.
    real(ccs_real), intent(in) :: side_length          !< The length of the side.

    type(ccs_mesh) :: mesh                             !< The resulting mesh.

    integer(ccs_int) :: start_global    ! The (global) starting index of a partition
    integer(ccs_int) :: end_global      ! The (global) last index of a partition
    integer(ccs_int) :: i               ! Loop counter
    integer(ccs_int) :: ii              ! Zero-indexed loop counter (simplifies some operations)
    integer(ccs_int) :: index_counter   ! Local index counter
    integer(ccs_int) :: face_counter    ! Cell-local face counter
    integer(ccs_int) :: comm_rank       ! The process ID within the parallel environment
    integer(ccs_int) :: comm_size       ! The size of the parallel environment

    integer(ccs_int) :: index_nb        ! The local index of a neighbour cell
    integer(ccs_int) :: global_index_nb ! The global index of a neighbour cell

    if(nx .eq. ny .and. ny .eq. nz) then !< @note Must be a cube (for now) @endnote

      select type (par_env)
      type is (parallel_environment_mpi)

        ! Set the global mesh parameters
        mesh%topo%global_num_cells = nx * ny * nz 
        mesh%geo%h = side_length / real(nx, ccs_real) !< @note Assumes cube @endnote

        ! Associate aliases to make code easier to read
        associate (nglobal => mesh%topo%global_num_cells, &
                  h => mesh%geo%h)

          ! Determine ownership range
          comm_rank = par_env%proc_id
          comm_size = par_env%num_procs
          start_global = global_start(nglobal, par_env%proc_id, par_env%num_procs)
          mesh%topo%local_num_cells = local_count(nglobal, par_env%proc_id, par_env%num_procs)
          end_global = start_global + (mesh%topo%local_num_cells - 1)

          ! Set max number of faces (constant, 6)
          mesh%topo%max_faces = 6_ccs_int
          
          ! Allocate mesh arrays
          allocate (mesh%topo%global_indices(mesh%topo%local_num_cells))
          allocate (mesh%topo%num_nb(mesh%topo%local_num_cells))
          allocate (mesh%topo%nb_indices(mesh%topo%max_faces, mesh%topo%local_num_cells))
          allocate (mesh%topo%face_indices(mesh%topo%max_faces, mesh%topo%local_num_cells))

          ! Initialise mesh arrays
          mesh%topo%num_nb(:) = mesh%topo%max_faces ! All cells have 6 neighbours (possibly ghost/boundary cells)

          ! Initalise neighbour indices
          mesh%topo%nb_indices(:,:) = 0_ccs_int

          ! First set the global index of local cells
          index_counter = 1_ccs_int
          do i = start_global, end_global
            mesh%topo%global_indices(index_counter) = i
            index_counter = index_counter + 1
          end do

          ! Assemble cells and faces
          ! XXX: Negative neighbour indices are used to indicate boundaries using the same numbering
          !      as cell-relative neighbour indexing, i.e.
          !        -1 = left boundary
          !        -2 = right boundary
          !        -3 = bottom boundary
          !        -4 = top boundary
          !        -5 = back_boundary
          !        -6 = front_boundary
          index_counter = 1_ccs_int ! Set local indexing starting from 1...n
          do i = start_global, end_global

            ii = i - 1_ccs_int

            ! Construct left (1) face/neighbour
            face_counter = left
            if (modulo(ii, nx) == 0_ccs_int) then
              index_nb = -left
              global_index_nb = -left
            else
              index_nb = index_counter - 1_ccs_int
              global_index_nb = i - 1_ccs_int
            end if
            call build_local_mesh_add_neighbour(index_counter, face_counter, index_nb, global_index_nb, mesh)

            ! Construct right (2) face/neighbour
            face_counter = right
            if (modulo(ii, nx) == (nx - 1_ccs_int)) then
              index_nb = -right
              global_index_nb = -right
            else
              index_nb = index_counter + 1_ccs_int
              global_index_nb = i + 1_ccs_int
            end if
            call build_local_mesh_add_neighbour(index_counter, face_counter, index_nb, global_index_nb, mesh)

            ! Construct bottom (3) face/neighbour
            face_counter = bottom
            if (modulo(ii/nx, ny) == 0_ccs_int) then
              index_nb = -bottom
              global_index_nb = -bottom
            else
              index_nb = index_counter - nx
              global_index_nb = i - nx
            end if
            call build_local_mesh_add_neighbour(index_counter, face_counter, index_nb, global_index_nb, mesh)

            ! Construct top (4) face/neighbour
            face_counter = top
            if (modulo(ii/nx, ny) == (ny - 1_ccs_int)) then
              index_nb = -top
              global_index_nb = -top
            else
              index_nb = index_counter + nx
              global_index_nb = i + nx
            end if
            call build_local_mesh_add_neighbour(index_counter, face_counter, index_nb, global_index_nb, mesh)

            ! Construct back (5) face/neighbour
            face_counter = back
            if ((ii / (nx * ny)) == nz - 1) then
              index_nb = -back
              global_index_nb = -back
            else
              index_nb = index_counter + nx * ny
              global_index_nb = i + nx * ny
            end if
            call build_local_mesh_add_neighbour(index_counter, face_counter, index_nb, global_index_nb, mesh)

            ! Construct front (6) face/neighbour
            face_counter = front
            if ((ii / (nx * ny)) == 0_ccs_int) then
              index_nb = -front
              global_index_nb = -front
            else
              index_nb = index_counter - nx * ny
              global_index_nb = i - nx * ny
            end if
            call build_local_mesh_add_neighbour(index_counter, face_counter, index_nb, global_index_nb, mesh)

            index_counter = index_counter + 1_ccs_int

          end do
        end associate

        ! print*,"Neighbour indices: ",mesh%neighbour_indices

        mesh%topo%total_num_cells = size(mesh%topo%global_indices)
        mesh%topo%halo_num_cells = mesh%topo%total_num_cells - mesh%topo%local_num_cells

        allocate (mesh%geo%x_p(ndim, mesh%topo%total_num_cells))
        allocate (mesh%geo%x_f(ndim, mesh%topo%max_faces, mesh%topo%local_num_cells))
        allocate (mesh%geo%volumes(mesh%topo%total_num_cells))
        allocate (mesh%geo%face_areas(mesh%topo%max_faces, mesh%topo%local_num_cells))
        allocate (mesh%geo%face_normals(ndim, mesh%topo%max_faces, mesh%topo%local_num_cells))

        mesh%geo%volumes(:) = mesh%geo%h**3 !< @note Mesh is cube @endnote
        mesh%geo%face_normals(:, :, :) = 0.0_ccs_real
        mesh%geo%x_p(:, :) = 0.0_ccs_real
        mesh%geo%x_f(:, :, :) = 0.0_ccs_real
        mesh%geo%face_areas(:, :) = mesh%geo%h 

        associate (h => mesh%geo%h)
          do i = 1_ccs_int, mesh%topo%total_num_cells
            ii = mesh%topo%global_indices(i)

            associate (x_p => mesh%geo%x_p(:, i))
              ! Set cell centre
              x_p(1) = (modulo(ii - 1, nx) + 0.5_ccs_real) * h
              x_p(2) = (modulo((ii - 1) / nx, ny) + 0.5_ccs_real) * h
              x_p(3) = (((ii - 1) / (nx * ny)) + 0.5_ccs_real) * h
            end associate
          end do

          do i = 1_ccs_int, mesh%topo%local_num_cells
            associate (x_p => mesh%geo%x_p(:, i), &
                       x_f => mesh%geo%x_f(:, :, i), &
                       normal => mesh%geo%face_normals(:, :, i))

              face_counter = left
              x_f(1, face_counter) = x_p(1) - 0.5_ccs_real * h
              x_f(2, face_counter) = x_p(2)
              normal(1, face_counter) = -1.0_ccs_real
              normal(2, face_counter) = 0.0_ccs_real
              normal(3, face_counter) = 0.0_ccs_real

              face_counter = right
              x_f(1, face_counter) = x_p(1) + 0.5_ccs_real * h
              x_f(2, face_counter) = x_p(2)
              normal(1, face_counter) = 1.0_ccs_real
              normal(2, face_counter) = 0.0_ccs_real
              normal(3, face_counter) = 0.0_ccs_real

              face_counter = bottom
              x_f(1, face_counter) = x_p(1)
              x_f(2, face_counter) = x_p(2) - 0.5_ccs_real * h
              normal(1, face_counter) = 0.0_ccs_real
              normal(2, face_counter) = -1.0_ccs_real
              normal(3, face_counter) = 0.0_ccs_real

              face_counter = top
              x_f(1, face_counter) = x_p(1)
              x_f(2, face_counter) = x_p(2) + 0.5_ccs_real * h
              normal(1, face_counter) = 0.0_ccs_real
              normal(2, face_counter) = 1.0_ccs_real
              normal(3, face_counter) = 0.0_ccs_real

              face_counter = back
              x_f(1, face_counter) = x_p(1)
              x_f(2, face_counter) = x_p(2) + 0.5_ccs_real * h
              normal(1, face_counter) = 0.0_ccs_real
              normal(2, face_counter) = 0.0_ccs_real
              normal(3, face_counter) = -1.0_ccs_real

              face_counter = front
              x_f(1, face_counter) = x_p(1)
              x_f(2, face_counter) = x_p(2) + 0.5_ccs_real * h
              normal(1, face_counter) = 0.0_ccs_real
              normal(2, face_counter) = 0.0_ccs_real
              normal(3, face_counter) = 1.0_ccs_real

            end associate
          end do
        end associate

        mesh%topo%num_faces = count_mesh_faces(mesh)

        call set_cell_face_indices(mesh)

      class default
        call error_abort("Unknown parallel environment type.")

      end select

    else
      print*,"Only supporting cubes for now - nx, ny and nz must be the same!"
    end if

  end function build_mesh

  !v Helper subroutine to add a neighbour to a cell's neighbour list.
  !
  !  Given a local and global index for a neighbour there are 3 possibilities:
  !
  !  1. the local and the neighbour is added immediately
  !  2. the global index is negative indicating it is a boundary and the "neighbour" is
  !     added immediately
  !  3. the index is not local:
  !     1. the global index is already in the off-process list (halos), the neighbour
  !        is added immediately
  !     2. this is a new halo cell, the list of global indices must be grown to
  !        accomodate before adding the neighbour.
  subroutine build_local_mesh_add_neighbour(index_p, index_p_nb, index_nb, global_index_nb, mesh)

    integer(ccs_int), intent(in) :: index_p !< the index of the cell whose neighbours we are assembling
    integer(ccs_int), intent(in) :: index_p_nb !< the cell-relative neighbour index
    integer(ccs_int), intent(in) :: index_nb !< the local index of the neighbour cell
    integer(ccs_int), intent(in) :: global_index_nb !< the global index of the neighbour cell
    type(ccs_mesh), intent(inout) :: mesh !< the mesh we are assembling neighbours on

    integer(ccs_int) :: ng  ! The current number of cells (total = local + halos)
    logical :: found        ! Indicates whether a halo cell was already present
    integer(ccs_int) :: i   ! Cell iteration counter

    if ((index_nb >= 1_ccs_int) .and. (index_nb <= mesh%topo%local_num_cells)) then
      ! Neighbour is local
      mesh%topo%nb_indices(index_p_nb, index_p) = index_nb
    else if (global_index_nb < 0_ccs_int) then
      ! Boundary "neighbour" - local index should also be -ve
      if (.not. (index_nb < 0_ccs_int)) then
        call error_abort("ERROR: boundary neighbours should have -ve indices.")
      end if
      mesh%topo%nb_indices(index_p_nb, index_p) = index_nb
    else
      ! Neighbour is in a halo

      ! First check if neighbour is already present in halo
      ng = size(mesh%topo%global_indices)
      found = .false.
      do i = mesh%topo%local_num_cells + 1, ng
        if (mesh%topo%global_indices(i) == global_index_nb) then
          found = .true.
          mesh%topo%nb_indices(index_p_nb, index_p) = i
          exit
        end if
      end do

      ! If neighbour was not present append to global index list (the end of the global index list
      ! becoming its local index).
      ! XXX: Note this currently copies into an n+1 temporary, reallocates and then copies back to
      !      the (extended) original array.
      if (.not. found) then
        if ((ng + 1) > mesh%topo%global_num_cells) then
          call error_abort("ERROR: Trying to create halo that exceeds global mesh size.")
        end if

        call append_to_arr(global_index_nb, mesh%topo%global_indices)
        ng = size(mesh%topo%global_indices)
        mesh%topo%nb_indices(index_p_nb, index_p) = ng
      end if
    end if

  end subroutine build_local_mesh_add_neighbour

  !v @note Docs needed.
  subroutine append_to_arr(i, arr)

    integer(ccs_int), intent(in) :: i
    integer(ccs_int), dimension(:), allocatable, intent(inout) :: arr ! XXX: Allocatable here be
                                                                      !      dragons! If this were
                                                                      !      intent(out) it would
                                                                      !      be deallocated on entry!
    integer(ccs_int) :: n
    integer(ccs_int), dimension(:), allocatable :: tmp

    n = size(arr)

    allocate (tmp(n + 1))

    tmp(1:n) = arr(1:n)

    n = n + 1
    tmp(n) = i

    deallocate (arr)
    allocate (arr(n))
    arr(:) = tmp(:)
    deallocate (tmp)

  end subroutine append_to_arr

  !v Count the number of faces in the mesh
  function count_mesh_faces(mesh) result(nfaces)

    ! Arguments
    type(ccs_mesh), intent(in) :: mesh !< the mesh

    ! Result
    integer(ccs_int) :: nfaces !< number of cell faces

    ! Local variables
    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb
    integer(ccs_int) :: global_index_p, index_p
    integer(ccs_int) :: j
    integer(ccs_int) :: nnb
    integer(ccs_int) :: n_faces_internal       ! Internal face count
    integer(ccs_int) :: nfaces_bnd       ! Boundary face count
    integer(ccs_int) :: nfaces_interface ! Process interface face count
    logical :: is_boundary
    logical :: is_local

    ! Initialise
    n_faces_internal = 0
    nfaces_bnd = 0
    nfaces_interface = 0

    ! Loop over cells
    do index_p = 1, mesh%topo%local_num_cells
      call set_cell_location(mesh, index_p, loc_p)
      call get_global_index(loc_p, global_index_p)
      call count_neighbours(loc_p, nnb)

      do j = 1, nnb
        call set_neighbour_location(loc_p, j, loc_nb)
        call get_boundary_status(loc_nb, is_boundary)

        if (.not. is_boundary) then
          call get_local_status(loc_nb, is_local)

          if (is_local) then
            ! Interior face
            n_faces_internal = n_faces_internal + 1
          else
            ! Process boundary face
            nfaces_interface = nfaces_interface + 1
          end if
        else
          ! Boundary face
          nfaces_bnd = nfaces_bnd + 1
        end if
      end do
    end do

    ! Interior faces will be counted twice
    nfaces = (n_faces_internal / 2) + nfaces_interface + nfaces_bnd

  end function count_mesh_faces

  !v @note Docs needed.
  subroutine set_cell_face_indices(mesh)

    ! Arguments
    type(ccs_mesh), intent(inout) :: mesh

    ! Local variables
    type(cell_locator) :: loc_p           ! Current cell
    type(neighbour_locator) :: loc_nb     ! Neighbour
    integer(ccs_int) :: index_nb, index_p
    integer(ccs_int) :: index_f
    integer(ccs_int) :: nnb
    integer(ccs_int) :: j
    integer(ccs_int) :: face_counter              ! Face index counter
    logical :: is_boundary

    face_counter = 0

    ! Loop over cells
    do index_p = 1, mesh%topo%local_num_cells
      call set_cell_location(mesh, index_p, loc_p)
      call count_neighbours(loc_p, nnb)

      do j = 1, nnb
        call set_neighbour_location(loc_p, j, loc_nb)
        call get_local_index(loc_nb, index_nb)
        call get_boundary_status(loc_nb, is_boundary)

        if (.not. is_boundary) then
          ! Cell with lowest local index assigns an index to the face
          if (index_p < index_nb) then
            face_counter = face_counter + 1
            call set_face_index(index_p, j, face_counter, mesh)
          else
            ! Find corresponding face in neighbour cell
            ! (To be improved, this seems inefficient!)
            index_f = get_neighbour_face_index(mesh, index_p, index_nb)
            call set_face_index(index_p, j, index_f, mesh)
          end if
        else
          face_counter = face_counter + 1
          call set_face_index(index_p, j, face_counter, mesh)
        end if
      end do  ! End loop over current cell's neighbours
    end do    ! End loop over local cells

  end subroutine set_cell_face_indices

  !v Computes the index of the face shared by the cells denoted by the specified
  !  local index and neighbouring index
  function get_neighbour_face_index(mesh, index_p, index_nb) result(index_f)
    type(ccs_mesh), intent(in) :: mesh       !< the mesh
    integer(ccs_int), intent(in) :: index_p  !< the current cell index
    integer(ccs_int), intent(in) :: index_nb !< the index of the neighbouring cell
    integer(ccs_int) :: index_f

    ! Local variables
    integer(ccs_int) :: k
    integer(ccs_int) :: nnb_nb
    type(cell_locator) :: loc_nb
    type(neighbour_locator) :: loc_nb_nb
    type(face_locator) :: loc_f
    integer(ccs_int) :: index_nb_nb

    call set_cell_location(mesh, index_nb, loc_nb)
    call count_neighbours(loc_nb, nnb_nb)
    do k = 1, nnb_nb
      call set_neighbour_location(loc_nb, k, loc_nb_nb)
      call get_local_index(loc_nb_nb, index_nb_nb)
      if (index_nb_nb == index_p) then
        call set_face_location(mesh, index_nb, k, loc_f)
        call get_local_index(loc_f, index_f)
        exit ! Exit the loop, as found shared face
      else if (k == nnb_nb) then
        call error_abort("ERROR: Failed to find face in owning cell.")
      end if
    end do
  end function get_neighbour_face_index

  integer function global_start(n, procid, nproc)

    integer(ccs_int), intent(in) :: n
    integer(ccs_int), intent(in) :: procid
    integer(ccs_int), intent(in) :: nproc

    ! Each PE gets an equal split of the problem with any remainder split equally between the lower
    ! PEs.
    global_start = procid * (n / nproc) + min(procid, modulo(n, nproc))

    ! Fortran indexing
    global_start = global_start + 1

  end function global_start

  integer function local_count(n, procid, nproc)

    integer(ccs_int), intent(in) :: n
    integer(ccs_int), intent(in) :: procid
    integer(ccs_int), intent(in) :: nproc

    if (procid < n) then
      local_count = global_start(n, procid, nproc)
      if (procid < (nproc - 1)) then
        local_count = global_start(n, procid + 1, nproc) - local_count
      else
        local_count = n - (local_count - 1)
      end if
    else
      local_count = 0
    end if

  end function local_count

end module mesh_utils
