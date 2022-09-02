module mesh_utils
#include "ccs_macros.inc"

  use constants, only: ndim, geoext, adiosconfig
  use utils, only: exit_print
  use kinds, only: ccs_int, ccs_long, ccs_real
  use types, only: ccs_mesh, topology, geometry, &
                   io_environment, io_process, &
                   face_locator, cell_locator, neighbour_locator
  use io, only: read_scalar, read_array, &
                configure_io, open_file, close_file, &
                initialise_io, cleanup_io
  use parallel, only: read_command_line_arguments
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
  !             4
  !       +----------+
  !       |          |
  !       |          |
  !     1 |          | 2
  !       |          |
  !       +----------+
  !             3
  !
  !  @endnote
  integer, parameter :: left = 1_ccs_int
  integer, parameter :: right = 2_ccs_int
  integer, parameter :: down = 3_ccs_int
  integer, parameter :: up = 4_ccs_int

  private
  public :: build_square_mesh
  public :: global_start
  public :: local_count
  public :: count_mesh_faces
  public :: read_mesh

contains

  !v Read mesh from file
  subroutine read_mesh(par_env, case_name, mesh)

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    character(len=:), allocatable :: case_name
    type(ccs_mesh), intent(inout) :: mesh                                   !< The mesh

    ! Local variables
    character(len=:), allocatable :: geo_file    ! Geo file name
    character(len=:), allocatable :: adios2_file ! ADIOS2 config file name

    class(io_environment), allocatable :: io_env
    class(io_process), allocatable :: geo_reader

    ! integer(ccs_long), dimension(1) :: sel_start
    ! integer(ccs_long), dimension(1) :: sel_count

    ! integer(ccs_long), dimension(2) :: sel2_start
    ! integer(ccs_long), dimension(2) :: sel2_count

    geo_file = case_name//geoext
    adios2_file = case_name//adiosconfig

    call initialise_io(par_env, adios2_file, io_env)
    call configure_io(io_env, "geo_reader", geo_reader)

    call open_file(geo_file, "read", geo_reader)

    call read_topology(par_env, geo_reader, mesh)
    call read_geometry(par_env, geo_reader, mesh)

    ! Close the file and ADIOS2 engine
    call close_file(geo_reader)

    ! Finalise the ADIOS2 IO environment
    call cleanup_io(io_env)

  end subroutine read_mesh

  !v Read the topology data from an input (HDF5) file
  ! This subroutine assumes the following names are used in the file:
  ! "ncel" - the total number of cells
  ! "nfac" - the total number of faces
  ! "maxfaces" - the maximum number of faces per cell
  ! "/face/cell1" and "/face/cell2" - the arrays the face edge data
  subroutine read_topology(par_env, geo_reader, mesh)

    use partitioning, only: compute_partitioner_input

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    class(io_process) :: geo_reader                                         !< The IO process for reading the file
    type(ccs_mesh), intent(inout) :: mesh                                   !< The mesh that will be read

    integer(ccs_int) :: i, j, k

    integer(ccs_long), dimension(1) :: sel_start
    integer(ccs_long), dimension(1) :: sel_count

    integer(ccs_long), dimension(2) :: sel2_start
    integer(ccs_long), dimension(2) :: sel2_count

    ! Read attribute "ncel" - the total number of cells
    call read_scalar(geo_reader, "ncel", mesh % topo % global_num_cells)
    ! Read attribute "nfac" - the total number of faces
    call read_scalar(geo_reader, "nfac", mesh % topo % global_num_faces)
    ! Read attribute "maxfaces" - the maximum number of faces per cell
    call read_scalar(geo_reader, "maxfaces", mesh % topo % max_faces)

    allocate (mesh % topo % face_cell1(mesh % topo % global_num_faces))
    allocate (mesh % topo % face_cell2(mesh % topo % global_num_faces))
    allocate (mesh % topo % global_face_indices(mesh % topo % max_faces, mesh % topo % global_num_cells))

    sel_start(1) = 0 ! Global index to start reading from
    sel_count(1) = mesh % topo % global_num_faces ! How many elements to read in total

    ! Read arrays face/cell1 and face/cell2
    call read_array(geo_reader, "/face/cell1", sel_start, sel_count, mesh % topo % face_cell1)
    call read_array(geo_reader, "/face/cell2", sel_start, sel_count, mesh % topo % face_cell2)

    sel2_start = 0
    sel2_count(1) = mesh % topo % max_faces! topo%global_num_cells
    sel2_count(2) = mesh % topo % global_num_cells

    call read_array(geo_reader, "/cell/cface", sel2_start, sel2_count, mesh % topo % global_face_indices)

    ! Create and populate the vtxdist array based on the total number of cells
    ! and the total number of ranks in the parallel environment
    allocate (mesh % topo % vtxdist(par_env % num_procs + 1)) ! vtxdist array is of size num_procs + 1 on all ranks

    mesh % topo % vtxdist(1) = 1                                                  ! First element is 1
    mesh % topo % vtxdist(par_env % num_procs + 1) = mesh % topo % global_num_cells + 1 ! Last element is total number of cells + 1

    ! Divide the total number of cells by the world size to
    ! compute the chunk sizes
    k = int(real(mesh % topo % global_num_cells) / par_env % num_procs)
    j = 1

    do i = 1, par_env % num_procs
      mesh % topo % vtxdist(i) = j
      j = j + k
    end do

    call compute_partitioner_input(par_env, mesh)

  end subroutine read_topology

  !v Read the geometry data from an input (HDF5) file
  subroutine read_geometry(par_env, geo_reader, mesh)

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    class(io_process) :: geo_reader                                         !< The IO process for reading the file
    type(ccs_mesh), intent(inout) :: mesh                                   !< The mesh%geometry that will be read

    integer(ccs_int) :: i, j, k, n, cell_count
    integer(ccs_int) :: start, end

    integer(ccs_long), dimension(1) :: vol_p_start
    integer(ccs_long), dimension(1) :: vol_p_count
    integer(ccs_long), dimension(2) :: x_p_start
    integer(ccs_long), dimension(2) :: x_p_count
    integer(ccs_long), dimension(2) :: f_xn_start
    integer(ccs_long), dimension(2) :: f_xn_count
    integer(ccs_long), dimension(1) :: f_a_start
    integer(ccs_long), dimension(1) :: f_a_count

    real(ccs_real), dimension(:, :), allocatable :: temp_x_f ! Temp array for face centres
    real(ccs_real), dimension(:, :), allocatable :: temp_n_f ! Temp array for face normals
    real(ccs_real), dimension(:), allocatable :: temp_a_f ! Temp array for face areas

    ! Read attribute "scalefactor"
    call read_scalar(geo_reader, "scalefactor", mesh % geo % scalefactor)

    ! Starting point for reading chunk of data
    vol_p_start = (/int(mesh % topo % vtxdist(par_env % proc_id + 1)) - 1/)
    ! How many data points will be read?
    vol_p_count = (/int(mesh % topo % vtxdist(par_env % proc_id + 2) - mesh % topo % vtxdist(par_env % proc_id + 1))/)

    ! Allocate memory for cell volumes array on each MPI rank
    allocate (mesh % geo % volumes(vol_p_count(1)))

    ! Read variable "/cell/vol"
    call read_array(geo_reader, "/cell/vol", vol_p_start, vol_p_count, mesh % geo % volumes)

    ! Starting point for reading chunk of data
    x_p_start = (/0, int(mesh % topo % vtxdist(par_env % proc_id + 1)) - 1/)
    ! How many data points will be read?
    x_p_count = (/ndim, int(mesh % topo % vtxdist(par_env % proc_id + 2) - mesh % topo % vtxdist(par_env % proc_id + 1))/)

    ! Allocate memory for cell centre coordinates array on each MPI rank
    allocate (mesh % geo % x_p(x_p_count(1), x_p_count(2)))

    ! Read variable "/cell/x"
    call read_array(geo_reader, "/cell/x", x_p_start, x_p_count, mesh % geo % x_p)

    ! Allocate temporary arrays for face centres, face normals and face areas
    allocate (temp_x_f(ndim, mesh % topo % global_num_faces))
    allocate (temp_n_f(ndim, mesh % topo % global_num_faces))
    allocate (temp_a_f(mesh % topo % global_num_faces))

    f_xn_start = 0
    f_xn_count(1) = ndim
    f_xn_count(2) = mesh % topo % global_num_faces

    ! Read variable "/face/x"
    call read_array(geo_reader, "/face/x", f_xn_start, f_xn_count, temp_x_f)
    ! Read variable "/face/n"
    call read_array(geo_reader, "/face/n", f_xn_start, f_xn_count, temp_n_f)

    f_a_start = 0
    f_a_count(1) = mesh % topo % global_num_faces

    ! Read variable "/face/area"
    call read_array(geo_reader, "/face/area", f_a_start, f_a_count, temp_a_f)

    ! Compute start and end points for local cells in global context
    start = mesh % topo % vtxdist(par_env % proc_id + 1)
    end = mesh % topo % vtxdist(par_env % proc_id + 2) - 1

    ! Allocate face centres, face normals and face areas arrays
    allocate (mesh % geo % x_f(ndim, mesh % topo % max_faces, mesh % topo % local_num_cells))
    allocate (mesh % geo % face_normals(ndim, mesh % topo % max_faces, mesh % topo % local_num_cells))
    allocate (mesh % geo % face_areas(mesh % topo % max_faces, mesh % topo % local_num_cells))

    cell_count = 1

    do k = start, end ! loop over cells owned by current process

      do j = 1, mesh % topo % max_faces ! loop over all faces for each cell
        n = mesh % topo % global_face_indices(j, k)

        do i = 1, ndim ! loop over dimensions
          ! Map from temp array to mesh
          mesh % geo % x_f(i, j, cell_count) = temp_x_f(i, n)
          mesh % geo % face_normals(i, j, cell_count) = temp_n_f(i, n)
        end do

        mesh % geo % face_areas(j, cell_count) = temp_a_f(n)

      end do
      cell_count = cell_count + 1

    end do

    ! Delete temp arrays
    deallocate (temp_x_f)
    deallocate (temp_n_f)
    deallocate (temp_a_f)

  end subroutine read_geometry

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
      mesh % topo % global_num_cells = cps**2
      mesh % geo % h = side_length / real(cps, ccs_real)

      ! Associate aliases to make code easier to read
      associate (nglobal => mesh % topo % global_num_cells, &
                 h => mesh % geo % h)

        ! Determine ownership range
        comm_rank = par_env % proc_id
        comm_size = par_env % num_procs
        start_global = global_start(nglobal, par_env % proc_id, par_env % num_procs)
        mesh % topo % local_num_cells = local_count(nglobal, par_env % proc_id, par_env % num_procs)
        end_global = start_global + (mesh % topo % local_num_cells - 1)

        ! Allocate mesh arrays
        allocate (mesh % topo % global_indices(mesh % topo % local_num_cells))
        allocate (mesh % topo % num_nb(mesh % topo % local_num_cells))
        allocate (mesh % topo % nb_indices(4, mesh % topo % local_num_cells))
        allocate (mesh % topo % face_indices(4, mesh % topo % local_num_cells))

        ! Initialise mesh arrays
        mesh % topo % num_nb(:) = 4_ccs_int ! All cells have 4 neighbours (possibly ghost/boundary cells)

        ! First set the global index of local cells
        index_counter = 1_ccs_int
        do i = start_global, end_global
          mesh % topo % global_indices(index_counter) = i
          index_counter = index_counter + 1
        end do

        ! Assemble cells and faces
        ! XXX: Negative neighbour indices are used to indicate boundaries using the same numbering
        !      as cell-relative neighbour indexing, i.e.
        !        -1 = left boundary
        !        -2 = right boundary
        !        -3 = down boundary
        !        -4 = up boundary
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

          ! Construct down (3) face/neighbour
          face_counter = down
          if ((ii / cps) == 0_ccs_int) then
            index_nb = -down
            global_index_nb = -down
          else
            index_nb = index_counter - cps
            global_index_nb = i - cps
          end if
          call build_local_mesh_add_neighbour(index_counter, face_counter, index_nb, global_index_nb, mesh)

          ! Construct up (4) face/neighbour
          face_counter = up
          if ((ii / cps) == (cps - 1_ccs_int)) then
            index_nb = -up
            global_index_nb = -up
          else
            index_nb = index_counter + cps
            global_index_nb = i + cps
          end if
          call build_local_mesh_add_neighbour(index_counter, face_counter, index_nb, global_index_nb, mesh)

          index_counter = index_counter + 1_ccs_int
        end do
      end associate

      mesh % topo % total_num_cells = size(mesh % topo % global_indices)
      mesh % topo % halo_num_cells = mesh % topo % total_num_cells - mesh % topo % local_num_cells

      allocate (mesh % geo % x_p(ndim, mesh % topo % total_num_cells))
      allocate (mesh % geo % x_f(ndim, 4, mesh % topo % local_num_cells)) !< @note Currently hardcoded as a 2D mesh. @endnote
      allocate (mesh % geo % volumes(mesh % topo % total_num_cells))
      allocate (mesh % geo % face_areas(4, mesh % topo % local_num_cells))
      allocate (mesh % geo % face_normals(ndim, 4, mesh % topo % local_num_cells)) ! Currently hardcoded as a 2D mesh.

      mesh % geo % volumes(:) = mesh % geo % h**2 !< @note Mesh is square and 2D @endnote
      mesh % geo % face_normals(:, :, :) = 0.0_ccs_real
      mesh % geo % x_p(:, :) = 0.0_ccs_real
      mesh % geo % x_f(:, :, :) = 0.0_ccs_real
      mesh % geo % face_areas(:, :) = mesh % geo % h  ! Mesh is square and 2D

      associate (h => mesh % geo % h)
        do i = 1_ccs_int, mesh % topo % total_num_cells
          ii = mesh % topo % global_indices(i)

          associate (x_p => mesh % geo % x_p(:, i))
            ! Set cell centre
            x_p(1) = (modulo(ii - 1, cps) + 0.5_ccs_real) * h
            x_p(2) = ((ii - 1) / cps + 0.5_ccs_real) * h
          end associate
        end do

        do i = 1_ccs_int, mesh % topo % local_num_cells
          associate (x_p => mesh % geo % x_p(:, i), &
                     x_f => mesh % geo % x_f(:, :, i), &
                     normal => mesh % geo % face_normals(:, :, i))

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

            face_counter = down
            x_f(1, face_counter) = x_p(1)
            x_f(2, face_counter) = x_p(2) - 0.5_ccs_real * h
            normal(1, face_counter) = 0.0_ccs_real
            normal(2, face_counter) = -1.0_ccs_real

            face_counter = up
            x_f(1, face_counter) = x_p(1)
            x_f(2, face_counter) = x_p(2) + 0.5_ccs_real * h
            normal(1, face_counter) = 0.0_ccs_real
            normal(2, face_counter) = 1.0_ccs_real
          end associate
        end do
      end associate

      mesh % topo % num_faces = count_mesh_faces(mesh)

      call set_cell_face_indices(mesh)

    class default
      call error_abort("Unknown parallel environment type.")

    end select
  end function build_square_mesh

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

    if ((index_nb >= 1_ccs_int) .and. (index_nb <= mesh % topo % local_num_cells)) then
      ! Neighbour is local
      mesh % topo % nb_indices(index_p_nb, index_p) = index_nb
    else if (global_index_nb < 0_ccs_int) then
      ! Boundary "neighbour" - local index should also be -ve
      if (.not. (index_nb < 0_ccs_int)) then
        call error_abort("ERROR: boundary neighbours should have -ve indices.")
      end if
      mesh % topo % nb_indices(index_p_nb, index_p) = index_nb
    else
      ! Neighbour is in a halo

      ! First check if neighbour is already present in halo
      ng = size(mesh % topo % global_indices)
      found = .false.
      do i = mesh % topo % local_num_cells + 1, ng
        if (mesh % topo % global_indices(i) == global_index_nb) then
          found = .true.
          mesh % topo % nb_indices(index_p_nb, index_p) = i
          exit
        end if
      end do

      ! If neighbour was not present append to global index list (the end of the global index list
      ! becoming its local index).
      ! XXX: Note this currently copies into an n+1 temporary, reallocates and then copies back to
      !      the (extended) original array.
      if (.not. found) then
        if ((ng + 1) > mesh % topo % global_num_cells) then
          call error_abort("ERROR: Trying to create halo that exceeds global mesh size.")
        end if

        call append_to_arr(global_index_nb, mesh % topo % global_indices)
        ng = size(mesh % topo % global_indices)
        mesh % topo % nb_indices(index_p_nb, index_p) = ng
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
    do index_p = 1, mesh % topo % local_num_cells
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
    do index_p = 1, mesh % topo % local_num_cells
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
