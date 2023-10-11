module mesh_utils
#include "ccs_macros.inc"

  use mpi

  use constants, only: ndim, geoext, adiosconfig
  use utils, only: exit_print, str, debug_print
  use kinds, only: ccs_int, ccs_long, ccs_real, ccs_err
  use types, only: ccs_mesh, topology, geometry, &
                   io_environment, io_process, &
                   face_locator, cell_locator, neighbour_locator, vert_locator, &
                   graph_connectivity
  use io, only: read_scalar, read_array, &
                write_scalar, write_array, &
                configure_io, open_file, close_file, &
                initialise_io, cleanup_io
  use parallel, only: read_command_line_arguments, create_shared_array, is_root, is_valid, create_shared_roots_comm, &
                      destroy_shared_array, sync
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi
  use meshing, only: get_global_index, get_natural_index, get_local_index, count_neighbours, &
                     create_cell_locator, create_neighbour_locator, create_face_locator, create_vert_locator, &
                     set_face_index, get_boundary_status, get_local_status, &
                     get_centre, set_centre, &
                     set_area, set_normal, get_face_normal, &
                     set_total_num_cells, get_total_num_cells, &
                     get_local_num_cells, set_local_num_cells, set_face_interpolation, &
                     get_global_num_cells, set_global_num_cells, &
                     set_halo_num_cells, &
                     get_global_num_faces, set_global_num_faces, &
                     get_num_faces, set_num_faces, &
                     get_max_faces, set_max_faces, &
                     get_vert_per_cell, set_vert_per_cell, &
                     get_vert_nb_per_cell, set_vert_nb_per_cell, &
                     get_global_num_vertices, set_global_num_vertices, &
                     set_face_interpolation, &
                     set_local_index, &
                     set_global_index, &
                     get_mesh_generated, set_mesh_generated
  use bc_constants
  use reordering, only: reorder_cells, compute_bandwidth

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

  integer, parameter :: front_bottom_left = 1_ccs_int
  integer, parameter :: front_bottom_right = 2_ccs_int
  integer, parameter :: front_top_right = 3_ccs_int
  integer, parameter :: front_top_left = 4_ccs_int
  integer, parameter :: back_bottom_left = 5_ccs_int
  integer, parameter :: back_bottom_right = 6_ccs_int
  integer, parameter :: back_top_right = 7_ccs_int
  integer, parameter :: back_top_left = 8_ccs_int
  integer, parameter :: front_left = 9_ccs_int
  integer, parameter :: front_right = 10_ccs_int
  integer, parameter :: front_bottom = 11_ccs_int
  integer, parameter :: front_top = 12_ccs_int
  integer, parameter :: middle_bottom_left = 13_ccs_int
  integer, parameter :: middle_bottom_right = 14_ccs_int
  integer, parameter :: middle_top_right = 15_ccs_int
  integer, parameter :: middle_top_left = 16_ccs_int
  integer, parameter :: back_left = 17_ccs_int
  integer, parameter :: back_right = 18_ccs_int
  integer, parameter :: back_bottom = 19_ccs_int
  integer, parameter :: back_top = 20_ccs_int

  private
  public :: build_square_mesh
  public :: build_square_topology
  public :: build_mesh
  public :: read_mesh
  public :: write_mesh
  public :: global_start
  public :: local_count
  public :: count_mesh_faces
  public :: set_cell_face_indices
  public :: compute_face_interpolation
  public :: partition_stride
  public :: print_topo
  public :: print_geo

contains

  !v Read mesh from file
  subroutine read_mesh(par_env, shared_env, case_name, mesh)

    use partitioning, only: compute_connectivity_get_local_cells, &
                            compute_partitioner_input

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The parallel environment
    character(len=:), allocatable :: case_name
    type(ccs_mesh), intent(inout) :: mesh                                   !< The mesh

    ! Local variables
    character(len=:), allocatable :: geo_file    ! Geo file name
    character(len=:), allocatable :: adios2_file ! ADIOS2 config file name

    class(parallel_environment), allocatable, target :: reader_env !< The reader parallel environment
    class(io_environment), allocatable :: io_env
    class(io_process), allocatable :: geo_reader

    call set_mesh_generated(.false., mesh)

    geo_file = case_name // "_mesh" // geoext
    adios2_file = case_name // adiosconfig

    call create_shared_roots_comm(par_env, shared_env, reader_env)

    call initialise_io(reader_env, adios2_file, io_env)
    call configure_io(io_env, "geo_reader", geo_reader)

    call open_file(geo_file, "read", geo_reader)

    call read_topology(par_env, shared_env, reader_env, geo_reader, mesh)

    call compute_partitioner_input(par_env, shared_env, mesh)

    call mesh_partition_reorder(par_env, shared_env, mesh)

    call read_geometry(shared_env, reader_env, geo_reader, mesh)

    ! Close the file and ADIOS2 engine
    call close_file(geo_reader)

    ! Finalise the ADIOS2 IO environment
    call cleanup_io(io_env)

    ! TODO: cleanup reader parallel environment

    call cleanup_topo(shared_env, mesh)

  end subroutine read_mesh

  !v Read the topology data from an input (HDF5) file
  ! This subroutine assumes the following names are used in the file:
  ! "ncel" - the total number of cells
  ! "nfac" - the total number of faces
  ! "maxfaces" - the maximum number of faces per cell
  ! "/face/cell1" and "/face/cell2" - the arrays the face edge data
  !
  ! This high-level interface zeroes the topology object contained by the mesh before calling the
  ! lower-level routine to read the topology object and building the vertext neighbours.
  subroutine read_topology(par_env, shared_env, reader_env, geo_reader, mesh)
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: reader_env !< The reader parallel environment
    class(io_process) :: geo_reader                                         !< The IO process for reading the file
    type(ccs_mesh), intent(inout) :: mesh                                   !< The mesh that will be read

    ! Zero scalar topology values to have known initial state
    call set_global_num_cells(0_ccs_int, mesh)
    call set_global_num_faces(0_ccs_int, mesh)
    call set_max_faces(0_ccs_int, mesh)
    call set_local_num_cells(0_ccs_int, mesh)
    call set_halo_num_cells(0_ccs_int, mesh)
    call set_total_num_cells(0_ccs_int, mesh)
    call set_num_faces(0_ccs_int, mesh)
    call set_vert_per_cell(0_ccs_int, mesh)
    call set_vert_nb_per_cell(0_ccs_int, mesh)
    call set_global_num_vertices(0_ccs_int, mesh)
    call read_topology_topo(par_env, shared_env, reader_env, geo_reader, mesh%topo)
    call build_vertex_neighbours(par_env, shared_env, mesh)

  end subroutine read_topology

  !v Read the topology data from an input (HDF5) file
  ! This subroutine assumes the following names are used in the file:
  ! "ncel" - the total number of cells
  ! "nfac" - the total number of faces
  ! "maxfaces" - the maximum number of faces per cell
  ! "/face/cell1" and "/face/cell2" - the arrays the face edge data
  !
  ! This lower-level subroutine works directly with the topology object, as indicated by the `_topo`
  ! suffix.
  subroutine read_topology_topo(par_env, shared_env, reader_env, geo_reader, topo)
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: reader_env !< The reader parallel environment
    class(io_process) :: geo_reader                                         !< The IO process for reading the file
    type(topology), intent(inout) :: topo                                   !< The mesh topology that will be read

    integer(ccs_int) :: i
    integer(ccs_int) :: num_bnd !< global number of boundary faces
    integer(ccs_int), dimension(:), allocatable :: bnd_rid, bnd_face
    integer(ccs_long), dimension(1) :: sel_start
    integer(ccs_long), dimension(1) :: sel_count

    integer(ccs_long), dimension(2) :: sel2_start
    integer(ccs_long), dimension(2) :: sel2_count

    integer(ccs_int) :: global_num_cells
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: global_num_faces
    integer(ccs_int) :: max_faces
    integer(ccs_int) :: vert_per_cell

    character(:), allocatable :: error_message

    integer(ccs_err) :: ierr
    integer :: shared_comm

    select type (shared_env)
    type is (parallel_environment_mpi)
      shared_comm = shared_env%comm
    class default
      shared_comm = -42
      call error_abort("Unsupported shared environment")
    end select

    ! Read attribute "ncel" - the total number of cells
    if (is_valid(reader_env)) then
      call read_scalar(geo_reader, "ncel", topo%global_num_cells)
    end if

    call MPI_Bcast(topo%global_num_cells, 1, MPI_INTEGER, 0, shared_comm, ierr)
    call get_global_num_cells(topo, global_num_cells)

    call read_topology_connectivity(shared_env, reader_env, geo_reader, &
                                    topo%global_num_faces, topo%max_faces, &
                                    topo%face_cell1, topo%face_cell1_window, topo%face_cell2, topo%face_cell2_window)
    call set_naive_distribution(par_env, global_num_cells, topo%graph_conn)

    ! XXX: <It should be possible to enter the partitioner here>

    ! Abort the execution if there a fewer global cells than MPI ranks
    if (topo%global_num_cells < par_env%num_procs) then
      error_message = "ERROR: Global number of cells < number of ranks. &
                      &Reduce the number of MPI ranks or use a bigger mesh."
      call error_abort(error_message)
    end if

    if (is_valid(reader_env)) then
      ! Read attribute "nvrt" - the total number of vertices
      call read_scalar(geo_reader, "nvrt", topo%global_num_vertices)

      ! Read attribute "nbnd" - the total number of boundary faces
      call read_scalar(geo_reader, "nbnd", num_bnd)
    end if
    call MPI_Bcast(topo%global_num_vertices, 1, MPI_INTEGER, 0, shared_comm, ierr)
    call MPI_Bcast(num_bnd, 1, MPI_INTEGER, 0, shared_comm, ierr)
    
    call get_max_faces(topo, max_faces)
    if (max_faces == 6) then ! if cell are hexes
      call set_vert_per_cell(8, topo) ! 8 vertices per cell
      call set_vert_nb_per_cell(20, topo)
    else
      call error_abort("Currently only supporting hex cells.")
    end if

    call get_global_num_faces(topo, global_num_faces)

    ! Read bnd data
    call create_shared_array(shared_env, global_num_faces, topo%bnd_rid, &
                             topo%bnd_rid_window)

    if (is_valid(reader_env)) then
      allocate (bnd_rid(num_bnd))
      allocate (bnd_face(num_bnd))

      sel_start(1) = 0 ! Global index to start reading from
      sel_count(1) = num_bnd ! How many elements to read in total
      call read_array(geo_reader, "/bnd/rid", sel_start, sel_count, bnd_rid)
      call read_array(geo_reader, "/bnd/face", sel_start, sel_count, bnd_face)

      ! make sure inside faces=0 and boundary faces are negative
      topo%bnd_rid(:) = 0_ccs_int
      topo%bnd_rid(bnd_face(:)) = -(bnd_rid(:) + 1_ccs_int)
    end if
    call sync(shared_env)

    ! Read global face and vertex indices
    call get_vert_per_cell(topo, vert_per_cell)
    call create_shared_array(shared_env, [max_faces, global_num_cells], topo%global_face_indices, &
                             topo%global_face_indices_window)
    call create_shared_array(shared_env, [vert_per_cell, global_num_cells], topo%global_vertex_indices, &
                             topo%global_vertex_indices_window)
                             
    if (is_valid(reader_env)) then
      sel2_start = 0
      sel2_count(1) = max_faces ! topo%global_num_cells
      sel2_count(2) = global_num_cells

      call read_array(geo_reader, "/cell/cface", sel2_start, sel2_count, topo%global_face_indices)

      sel2_count(1) = vert_per_cell

      call read_array(geo_reader, "/cell/vertices", sel2_start, sel2_count, topo%global_vertex_indices)
    end if

    associate (irank => par_env%proc_id)
      local_num_cells = int(topo%graph_conn%vtxdist(irank + 2) - topo%graph_conn%vtxdist(irank + 1), ccs_int)
      call set_local_num_cells(local_num_cells, topo)
      call set_total_num_cells(local_num_cells, topo)

      allocate (topo%global_indices(local_num_cells))
      do i = 1, topo%local_num_cells
        topo%global_indices(i) = int(topo%graph_conn%vtxdist(irank + 1), ccs_int) + (i - 1)
      end do

      allocate (topo%num_nb(local_num_cells))
      topo%num_nb(:) = max_faces
    end associate

  end subroutine read_topology_topo

  subroutine read_topology_connectivity(shared_env, reader_env, geo_reader, &
                                        global_num_faces, max_faces, &
                                        face_cell1, face_cell1_window, face_cell2, face_cell2_window)

    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: reader_env !< The reader parallel environment
    class(io_process) :: geo_reader                                            !< The IO process for reading the file
    integer(ccs_int), intent(out) :: global_num_faces
    integer(ccs_int), intent(out) :: max_faces
    integer(ccs_int), dimension(:), pointer, intent(out) :: face_cell1
    integer, intent(out) :: face_cell1_window
    integer(ccs_int), dimension(:), pointer, intent(out) :: face_cell2
    integer, intent(out) :: face_cell2_window

    integer(ccs_long), dimension(1) :: sel_start
    integer(ccs_long), dimension(1) :: sel_count

    integer :: shared_comm
    integer(ccs_err) :: ierr

    select type (shared_env)
    type is (parallel_environment_mpi)
      shared_comm = shared_env%comm
    class default
      shared_comm = -42
      call error_abort("Unsupported shared environment")
    end select

    if (is_valid(reader_env)) then
      ! Read attribute "nfac" - the total number of faces
      call read_scalar(geo_reader, "nfac", global_num_faces)
      ! Read attribute "maxfaces" - the maximum number of faces per cell
      call read_scalar(geo_reader, "maxfaces", max_faces)
    end if
    call MPI_Bcast(global_num_faces, 1, MPI_INTEGER, 0, shared_comm, ierr)
    call MPI_Bcast(max_faces, 1, MPI_INTEGER, 0, shared_comm, ierr)

    ! Read arrays face/cell1 and face/cell2
    call create_shared_array(shared_env, global_num_faces, face_cell1, &
                             face_cell1_window)
    call create_shared_array(shared_env, global_num_faces, face_cell2, &
                             face_cell2_window)

    if (is_valid(reader_env)) then
      sel_start = 0                ! Global index to start reading from
      sel_count = global_num_faces ! How many elements to read in total
      call read_array(geo_reader, "/face/cell1", sel_start, sel_count, face_cell1)
      call read_array(geo_reader, "/face/cell2", sel_start, sel_count, face_cell2)
    end if
    call sync(shared_env)

  end subroutine read_topology_connectivity

  !v Build the vertex neighbours from the cell-vertex connectivity.
  !
  !  @note@ This will be quadratic - do we REALLY need it?
  !  @note@ This won't find boundary vertex "neighbours" hopefully they aren't needed...
  subroutine build_vertex_neighbours(par_env, shared_env, mesh)

    use case_config, only: vertex_neighbours

    class(parallel_environment), intent(in) :: par_env
    class(parallel_environment), intent(in) :: shared_env
    type(ccs_mesh), intent(inout) :: mesh

    integer(ccs_int), dimension(:), pointer :: global_num_vert_nb           !< The local number of vertex neighbours per cell
    integer :: global_num_vert_nb_window                                    !< Associated shared window
    integer(ccs_int) :: global_num_cells
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: vert_nb_per_cell
    integer(ccs_int) :: vert_per_cell

    integer(ccs_int) :: i, j, k
    integer(ccs_int) :: global_vert_index

    integer(ccs_int) :: idx_vnb
    integer(ccs_int) :: total_num_cells

    integer(ccs_int) :: vctr

    associate (foo => shared_env)
    end associate
    if (vertex_neighbours .eqv. .true.) then

      if (par_env%proc_id == par_env%root) then
        print *, "Building vertex neighbours, this may take a while..."
      end if

      call get_global_num_cells(mesh, global_num_cells)
      call get_vert_nb_per_cell(mesh, vert_nb_per_cell)
      call get_vert_per_cell(mesh, vert_per_cell)

      call create_shared_array(shared_env, [vert_nb_per_cell, global_num_cells], mesh%topo%global_vert_nb_indices, mesh%topo%global_vert_nb_indices_window)
      call create_shared_array(shared_env, global_num_cells, global_num_vert_nb, global_num_vert_nb_window) ! XXX: Will have to shrink this later...

      if (is_root(shared_env)) then
        mesh%topo%global_vert_nb_indices(:, :) = 0 ! Not an internal neighbour, not a boundary - will this work?

        ! Yes, quadratic...
        do i = 1, global_num_cells
          global_num_vert_nb(i) = 0
          do j = 1, global_num_cells
            if (j /= i) then
              if (any(mesh%topo%global_face_indices(:, i) == j)) then
                ! Face neighbour, ignore
                continue
              else
                do k = 1, vert_per_cell
                  global_vert_index = mesh%topo%global_vertex_indices(k, i)
                  if (any(mesh%topo%global_vertex_indices(:, j) == global_vert_index)) then
                    global_num_vert_nb(i) = global_num_vert_nb(i) + 1
                    mesh%topo%global_vert_nb_indices(global_num_vert_nb(i), i) = j
                    exit
                  end if
                end do
              end if
            end if

            ! Found all the vertex neighbours we're expecting
            if (global_num_vert_nb(i) == vert_per_cell) then
              exit
            end if
          end do
        end do

      end if
      call sync(shared_env)

      ! Get local data
      call get_local_num_cells(mesh, local_num_cells)
      allocate (mesh%topo%num_vert_nb(local_num_cells))
      allocate (mesh%topo%vert_nb_indices(vert_nb_per_cell, local_num_cells))
      do i = 1, local_num_cells
        vctr = 1
        associate (idxg => mesh%topo%global_indices(i))
          mesh%topo%num_vert_nb(i) = global_num_vert_nb(idxg)
          mesh%topo%vert_nb_indices(:, i) = mesh%topo%global_vert_nb_indices(:, idxg)
          do j = 1, mesh%topo%num_vert_nb(i)
            associate (idxg_vnb => mesh%topo%vert_nb_indices(j, i))
              if (idxg_vnb > 0) then
                idx_vnb = findloc(mesh%topo%global_indices, idxg_vnb, dim=1)

                if (idx_vnb > 0) then
                  call build_local_mesh_add_neighbour(i, vctr, idx_vnb, idxg_vnb, mesh, .true.)
                else
                  call get_total_num_cells(mesh, total_num_cells)
                  call build_local_mesh_add_neighbour(i, vctr, total_num_cells + 1, idxg_vnb, mesh, .true.)
                end if

                vctr = vctr + 1
              end if
            end associate
          end do
        end associate
      end do

      call get_total_num_cells(mesh, total_num_cells)
      if (any(mesh%topo%vert_nb_indices > total_num_cells)) then
        call error_abort("ERROR: Vertex neighbour index outside total number of cells I can see")
      end if

      call destroy_shared_array(shared_env, mesh%topo%global_vert_nb_indices, mesh%topo%global_vert_nb_indices_window)
      call destroy_shared_array(shared_env, global_num_vert_nb, global_num_vert_nb_window)

    else
      if (par_env%proc_id == par_env%root) then
        print *, "Not building vertex neighbours"
      end if

    end if

  end subroutine build_vertex_neighbours

  !v Read the geometry data from an input (HDF5) file
  subroutine read_geometry(shared_env, reader_env, geo_reader, mesh)

    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: reader_env !< The reader parallel environment
    class(io_process) :: geo_reader                                         !< The IO process for reading the file
    type(ccs_mesh), intent(inout) :: mesh                                   !< The mesh%geometry that will be read

    integer(ccs_int) :: i, j, n, global_icell, local_icell
    integer(ccs_int) :: vert_per_cell

    integer(ccs_long), dimension(1) :: vol_p_start
    integer(ccs_long), dimension(1) :: vol_p_count
    integer(ccs_long), dimension(2) :: x_p_start
    integer(ccs_long), dimension(2) :: x_p_count
    integer(ccs_long), dimension(2) :: f_xn_start
    integer(ccs_long), dimension(2) :: f_xn_count
    integer(ccs_long), dimension(1) :: f_a_start
    integer(ccs_long), dimension(1) :: f_a_count

    real(ccs_real), dimension(:), pointer :: temp_vol_c ! Temp array for cell volumes
    real(ccs_real), dimension(:, :), pointer :: temp_x_p ! Temp array for cell centres
    real(ccs_real), dimension(:, :), pointer :: temp_x_f ! Temp array for face centres
    real(ccs_real), dimension(:, :), pointer :: temp_n_f ! Temp array for face normals
    real(ccs_real), dimension(:, :), pointer :: temp_x_v ! Temp array for vertex coordinates
    real(ccs_real), dimension(:), pointer :: temp_a_f ! Temp array for face areas

    real(ccs_real), dimension(3) :: face_normal, x_p, x_f
    integer(ccs_int) :: local_num_cells, index_p, nnb
    type(cell_locator) :: loc_p ! Face locator object
    integer(ccs_int) :: total_num_cells
    integer(ccs_int) :: global_num_cells
    integer(ccs_int) :: global_num_faces
    integer(ccs_int) :: global_num_vertices
    integer(ccs_int) :: max_faces
    type(face_locator) :: loc_f ! Face locator object
    type(vert_locator) :: loc_v ! Vertex locator object

    integer(ccs_err) :: ierr
    integer :: shared_comm

    integer :: temp_a_f_window, temp_n_f_window, temp_window, temp_x_f_window, temp_x_v_window

    select type (shared_env)
    type is (parallel_environment_mpi)
      shared_comm = shared_env%comm
    class default
      shared_comm = -42
      call error_abort("Unsupported shared environment")
    end select

    call get_max_faces(mesh, max_faces)
    if (max_faces == 6) then ! if cell are hexes
      call set_vert_per_cell(8, mesh) ! 8 vertices per cell
    else
      call error_abort("Currently only supporting hex cells.")
    end if

    call get_vert_per_cell(mesh, vert_per_cell)

    ! Read attribute "scalefactor"
    if (is_valid(reader_env)) then
      call read_scalar(geo_reader, "scalefactor", mesh%geo%scalefactor)
    end if
    call MPI_Bcast(mesh%geo%scalefactor, 1, MPI_DOUBLE_PRECISION, 0, shared_comm, ierr)

    ! Starting point for reading chunk of data
    vol_p_start = 0

    ! How many data points will be read?
    call get_global_num_cells(mesh, global_num_cells)
    vol_p_count = global_num_cells

    ! Allocate memory for cell volumes array on each MPI rank
    call get_total_num_cells(mesh, total_num_cells)
    allocate (mesh%geo%volumes(total_num_cells))

    ! Read variable "/cell/vol"
    call create_shared_array(shared_env, global_num_cells, temp_vol_c, &
                             temp_window)
    if (is_valid(reader_env)) then
      call read_array(geo_reader, "/cell/vol", vol_p_start, vol_p_count, temp_vol_c)
    end if
    call sync(shared_env)
    mesh%geo%volumes(:) = temp_vol_c(mesh%topo%natural_indices(:))
    call sync(shared_env)
    call destroy_shared_array(shared_env, temp_vol_c, temp_window)

    ! Starting point for reading chunk of data
    x_p_start = [0, 0]

    ! How many data points will be read?
    x_p_count = [ndim, global_num_cells]

    ! Allocate memory for cell centre coordinates array on each MPI rank
    allocate (mesh%geo%x_p(ndim, total_num_cells))

    ! Read variable "/cell/x"
    call create_shared_array(shared_env, [ndim, global_num_cells], temp_x_p, &
                             temp_window)
    if (is_valid(reader_env)) then
      call read_array(geo_reader, "/cell/x", x_p_start, x_p_count, temp_x_p)
    end if
    call sync(shared_env)
    mesh%geo%x_p(:, :) = temp_x_p(:, mesh%topo%natural_indices(:))
    call sync(shared_env)
    call destroy_shared_array(shared_env, temp_x_p, temp_window)

    ! Allocate temporary arrays for face centres, face normals, face areas and vertex coords
    call get_global_num_faces(mesh, global_num_faces)

    call create_shared_array(shared_env, [ndim, global_num_faces], temp_x_f, &
                             temp_x_f_window)
    call create_shared_array(shared_env, [ndim, global_num_faces], temp_n_f, &
                             temp_n_f_window)
    call create_shared_array(shared_env, global_num_faces, temp_a_f, &
                             temp_a_f_window)

    f_xn_start = 0
    f_xn_count(1) = ndim
    f_xn_count(2) = global_num_faces

    if (is_valid(reader_env)) then
      ! Read variable "/face/x"
      call read_array(geo_reader, "/face/x", f_xn_start, f_xn_count, temp_x_f)
      ! Read variable "/face/n"
      call read_array(geo_reader, "/face/n", f_xn_start, f_xn_count, temp_n_f)
    end if

    f_a_start = 0
    f_a_count(1) = global_num_faces

    if (is_valid(reader_env)) then
      ! Read variable "/face/area"
      call read_array(geo_reader, "/face/area", f_a_start, f_a_count, temp_a_f)
    end if
    call sync(shared_env)

    ! Read variable "/vert"
    call get_global_num_vertices(mesh, global_num_vertices)
    call create_shared_array(shared_env, [ndim, global_num_vertices], temp_x_v, &
                             temp_x_v_window)
    f_xn_count(1) = ndim
    f_xn_count(2) = global_num_vertices
    if (is_valid(reader_env)) then
      call read_array(geo_reader, "/vert", f_xn_start, f_xn_count, temp_x_v)
    end if
    call sync(shared_env)

    ! Allocate arrays for face centres, face normals, face areas arrand vertex coordinates
    call get_local_num_cells(mesh, local_num_cells)
    allocate (mesh%geo%x_f(ndim, max_faces, local_num_cells))
    allocate (mesh%geo%face_normals(ndim, max_faces, local_num_cells))
    allocate (mesh%geo%face_areas(max_faces, local_num_cells))
    allocate (mesh%geo%vert_coords(ndim, vert_per_cell, local_num_cells))

    ! Procs fill local data
    call sync(shared_env)
    do local_icell = 1, local_num_cells ! loop over cells owned by current process
      call create_cell_locator(mesh, local_icell, loc_p)
      call get_natural_index(loc_p, global_icell)

      do j = 1, max_faces ! loop over all faces for each cell
        call create_face_locator(mesh, local_icell, j, loc_f)

        n = mesh%topo%global_face_indices(j, global_icell)
        call set_centre(loc_f, temp_x_f(:, n))
        do i = 1, ndim ! loop over dimensions
          ! Map from temp array to mesh for face centres and face normals
          mesh%geo%face_normals(i, j, local_icell) = temp_n_f(i, n)
        end do

        ! Map from temp array to mesh for face areas
        call set_area(temp_a_f(n), loc_f)
      end do

      do j = 1, vert_per_cell ! loop over all vertices for each cell
        call create_vert_locator(mesh, local_icell, j, loc_v)

        n = mesh%topo%loc_global_vertex_indices(j, local_icell)
        call set_centre(loc_v, temp_x_v(:, n))
      end do

    end do

    ! Correct normal orientations and norms
    do index_p = 1, local_num_cells ! loop over cells owned by current process

      call create_cell_locator(mesh, index_p, loc_p)
      call get_centre(loc_p, x_p)
      call count_neighbours(loc_p, nnb)

      do j = 1, nnb ! loop over all faces for each cell

        call create_face_locator(mesh, index_p, j, loc_f)
        call get_face_normal(loc_f, face_normal)
        call get_centre(loc_f, x_f)

        if (dot_product(face_normal(:), x_f - x_p) < 0.0_ccs_real) then
          face_normal = -face_normal
        end if

        ! Normalise face normals too
        call set_normal(loc_f, face_normal / norm2(face_normal))

      end do
    end do

    ! Delete temp arrays
    call sync(shared_env)
    call destroy_shared_array(shared_env, temp_x_f, temp_x_f_window)
    call destroy_shared_array(shared_env, temp_x_v, temp_x_v_window)
    call destroy_shared_array(shared_env, temp_n_f, temp_n_f_window)
    call destroy_shared_array(shared_env, temp_a_f, temp_a_f_window)

    call compute_face_interpolation(mesh)

  end subroutine read_geometry

  !v Write mesh to file
  subroutine write_mesh(par_env, case_name, mesh)

    ! Arguments
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    character(len=:), allocatable, intent(in) :: case_name                  !< The case name
    type(ccs_mesh), intent(inout) :: mesh                                   !< The mesh

    ! Local variables
    character(len=:), allocatable :: geo_file    ! Geo file name
    character(len=:), allocatable :: adios2_file ! ADIOS2 config file name

    class(io_environment), allocatable :: io_env
    class(io_process), allocatable :: geo_writer

    logical :: is_generated

    call get_mesh_generated(mesh, is_generated)

    if (.not. is_generated) then
      ! Mesh was read, no need to write again
      return
    end if

    ! Set ADIOS2 config file name
    adios2_file = case_name // adiosconfig

    ! Set geo file name
    geo_file = case_name // geoext

    ! Open geo file for writing
    call initialise_io(par_env, adios2_file, io_env)
    call configure_io(io_env, "geo_writer", geo_writer)
    call open_file(geo_file, "write", geo_writer)

    ! Write mesh
    call write_topology(par_env, geo_writer, mesh)
    call write_geometry(par_env, geo_writer, mesh)

    ! Close the file and ADIOS2 engine
    call close_file(geo_writer)

    ! Finalise the ADIOS2 environment
    call cleanup_io(io_env)

  end subroutine write_mesh

  !v Write the mesh topology data to file
  subroutine write_topology(par_env, geo_writer, mesh)

    use mpi

    ! Arguments
    class(parallel_environment), intent(in) :: par_env               !< The parallel environment
    class(io_process), allocatable, target, intent(in) :: geo_writer !< The IO process for writing the mesh ("geo") file
    type(ccs_mesh), intent(in) :: mesh                               !< The mesh

    ! Local variables
    integer(ccs_long), dimension(2) :: sel2_shape
    integer(ccs_long), dimension(2) :: sel2_start
    integer(ccs_long), dimension(2) :: sel2_count

    integer(ccs_int) :: global_num_cells
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: vert_per_cell

    type(cell_locator) :: loc_p
    integer(ccs_int) :: index_global

    integer(ccs_int), dimension(:), allocatable :: natural_vertices_1d
    integer(ccs_int), dimension(:, :), allocatable :: natural_vertices_2d
    integer(ccs_err) ierr

    integer(ccs_int) :: i, j, idx

    call get_local_num_cells(mesh, local_num_cells)
    call get_global_num_cells(mesh, global_num_cells)
    call get_vert_per_cell(mesh, vert_per_cell)

    if (vert_per_cell == 0) then
      call error_abort("Number of vertices per cell unset.")
    end if

    call create_cell_locator(mesh, 1, loc_p)
    call get_global_index(loc_p, index_global)

    ! Write cell vertices
    sel2_shape(1) = vert_per_cell
    sel2_shape(2) = global_num_cells
    sel2_start(1) = 0
    sel2_start(2) = index_global - 1
    sel2_count(1) = vert_per_cell
    sel2_count(2) = local_num_cells

    ! Get global vertex indices in cell natural order
    allocate (natural_vertices_1d(vert_per_cell * global_num_cells))
    natural_vertices_1d(:) = 0
    do i = 1, local_num_cells
      idx = vert_per_cell * (mesh%topo%natural_indices(i) - 1)
      do j = 1, vert_per_cell
        natural_vertices_1d(idx + j) = mesh%topo%loc_global_vertex_indices(j, i)
      end do
    end do
    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_Allreduce(MPI_IN_PLACE, natural_vertices_1d, size(natural_vertices_1d), &
                         MPI_INTEGER, MPI_SUM, par_env%comm, ierr)
    class default
      call error_abort("Unknown parallel environment")
    end select

    allocate (natural_vertices_2d(vert_per_cell, local_num_cells))
    do i = 1, local_num_cells
      idx = vert_per_cell * (mesh%topo%global_indices(i) - 1)
      do j = 1, vert_per_cell
        natural_vertices_2d(j, i) = natural_vertices_1d(idx + j)
      end do
    end do

    deallocate (natural_vertices_1d)

    call write_array(geo_writer, "/cell/vertices", sel2_shape, sel2_start, sel2_count, &
                     natural_vertices_2d)

    deallocate (natural_vertices_2d)

  end subroutine write_topology

  !v Write the mesh geometry data to file
  subroutine write_geometry(par_env, geo_writer, mesh)

    ! Arguments
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    class(io_process), allocatable, target, intent(in) :: geo_writer        !< The IO process for writing the mesh ("geo") file
    type(ccs_mesh), intent(inout) :: mesh                                   !< The mesh

    ! Local variables
    integer(ccs_long), dimension(2) :: sel2_shape
    integer(ccs_long), dimension(2) :: sel2_start
    integer(ccs_long), dimension(2) :: sel2_count

    integer(ccs_int) :: i
    integer(ccs_int) :: verts_per_side
    integer(ccs_int) :: vert_per_cell
    integer(ccs_int) :: global_num_vertices

    real(ccs_real), dimension(:, :), allocatable :: vert_coords_tmp

    call get_vert_per_cell(mesh, vert_per_cell)
    call get_global_num_vertices(mesh, global_num_vertices)

    if (vert_per_cell == 0) then
      call error_abort("Number of vertices per cell unset.")
    end if

    ! Root process calculates all (global) vertex coords and stores temporarily
    if (par_env%proc_id == par_env%root) then
      allocate (vert_coords_tmp(ndim, global_num_vertices))

      vert_coords_tmp = 0.0_ccs_real

      if (vert_per_cell == 4) then
        verts_per_side = nint(global_num_vertices**(1./2))
      else
        verts_per_side = nint(global_num_vertices**(1./3))
      end if

      do i = 1, global_num_vertices
        vert_coords_tmp(1, i) = modulo(i - 1, verts_per_side) * mesh%geo%h
        vert_coords_tmp(2, i) = modulo((i - 1) / verts_per_side, verts_per_side) * mesh%geo%h
        vert_coords_tmp(3, i) = ((i - 1) / (verts_per_side * verts_per_side)) * mesh%geo%h
      end do
    end if

    sel2_shape(1) = ndim
    sel2_shape(2) = global_num_vertices
    sel2_start(1) = 0
    sel2_start(2) = 0
    if (par_env%proc_id == par_env%root) then
      sel2_count(1) = ndim
      sel2_count(2) = global_num_vertices
    else
      sel2_count(1) = 0
      sel2_count(2) = 0
    end if

    call write_array(geo_writer, "/vert", sel2_shape, sel2_start, sel2_count, vert_coords_tmp)

    if (par_env%proc_id == par_env%root) then
      deallocate (vert_coords_tmp)
    end if

  end subroutine write_geometry

  !v Utility constructor to build a 2D mesh with hex cells.
  !
  !  Builds a Cartesian grid of nx*ny cells.
  function build_square_mesh(par_env, shared_env, cps, side_length) result(mesh)

    use partitioning, only: compute_partitioner_input

    class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared memory environment
    integer(ccs_int), intent(in) :: cps                !< Number of cells per side of the mesh.
    real(ccs_real), intent(in) :: side_length          !< The length of the side.

    type(ccs_mesh) :: mesh                             !< The resulting mesh.

    character(:), allocatable :: error_message

    if (cps * cps < par_env%num_procs) then
      error_message = "ERROR: Global number of cells < number of ranks. &
                      &Increase the mesh size or reduce the number of MPI ranks."
      call error_abort(error_message)
    end if

    call set_mesh_generated(.true., mesh)

    call build_square_topology(par_env, shared_env, cps, mesh)

    call compute_partitioner_input(par_env, shared_env, mesh)

    call mesh_partition_reorder(par_env, shared_env, mesh)

    call build_square_geometry(par_env, cps, side_length, mesh)

    call cleanup_topo(shared_env, mesh)

  end function build_square_mesh

  !v Utility constructor to build a square mesh.
  !
  !  Builds a Cartesian grid of NxN cells on the domain LxL.
  subroutine build_square_topology(par_env, shared_env, cps, mesh)

    class(parallel_environment), intent(in) :: par_env    !< The parallel environment to construct the mesh.
    class(parallel_environment), intent(in) :: shared_env !< The shared memory environment
    integer(ccs_int), intent(in) :: cps                !< Number of cells per side of the mesh.

    type(ccs_mesh), intent(inout) :: mesh                             !< The resulting mesh.

    integer(ccs_int) :: start_global    ! The (global) starting index of a partition
    integer(ccs_int) :: end_global      ! The (global) last index of a partition
    integer(ccs_int) :: i               ! Loop counter
    integer(ccs_int) :: ii              ! Zero-indexed loop counter (simplifies some operations)
    integer(ccs_int) :: index_counter   ! Local index counter
    integer(ccs_int) :: face_counter    ! Cell-local face counter
    integer(ccs_int) :: vertex_counter    ! Cell-local face counter
    integer(ccs_int) :: face_index_counter    ! global face counter

    integer(ccs_int) :: global_index_nb ! The global index of a neighbour cell

    integer(ccs_int) :: nglobal          ! The global number of cells
    integer(ccs_int) :: local_num_cells  ! The local number of cells
    integer(ccs_int) :: total_num_cells  ! The total number of cells
    integer(ccs_int) :: global_num_faces ! The global number of faces
    integer(ccs_int) :: max_faces        ! The maximum number of faces per cell
    integer(ccs_int) :: vert_per_cell    ! The number of vertices per cell
    integer(ccs_int) :: vert_nb_per_cell ! The number of neighbours via vertices per cell
    integer(ccs_int) :: global_num_vertices
    logical :: set_vert_nb               ! Flag for setting vertex neighbour
    integer(ccs_int), dimension(2) :: nb_direction  ! Array indicating direction of neighbour

    type(face_locator) :: loc_f

    integer(ccs_int), dimension(2) :: length

    nglobal = cps**2 ! The global cell count
    call build_square_topology_connectivity(shared_env, &
                                            cps, &
                                            global_num_faces, max_faces, &
                                            mesh%topo%face_cell1, mesh%topo%face_cell1_window, &
                                            mesh%topo%face_cell2, mesh%topo%face_cell2_window)
    call set_naive_distribution(par_env, nglobal, mesh%topo%graph_conn)

    ! XXX: <It should be possible to enter the partitioner here>

    select type (par_env)
    type is (parallel_environment_mpi)

      select type (shared_env)
      type is (parallel_environment_mpi)

        ! Set the global mesh parameters
        call set_global_num_cells(nglobal, mesh)
        call set_global_num_vertices((cps + 1)**2, mesh)

        call set_global_num_faces(global_num_faces, mesh)
        call set_max_faces(max_faces, mesh)

        ! Just to make sure we are working with the same numbers as the mesh object.
        call get_global_num_cells(mesh, nglobal)
        call get_global_num_vertices(mesh, global_num_vertices)
        call get_global_num_faces(mesh, global_num_faces)
        call get_max_faces(mesh, max_faces)

        ! Associate aliases to make code easier to read
        associate (h => mesh%geo%h)

          ! Determine ownership range
          start_global = global_start(nglobal, par_env%proc_id, par_env%num_procs)
          local_num_cells = local_count(nglobal, par_env%proc_id, par_env%num_procs)
          call set_local_num_cells(local_num_cells, mesh)
          call get_local_num_cells(mesh, local_num_cells) ! Ensure using correct value

          ! Abort the execution if any rank has 0 local cells
          if (local_num_cells <= 0) then
            call error_abort("ERROR: Zero local cells found.")
          end if

          call set_total_num_cells(local_num_cells, mesh) ! Set initial value
          end_global = start_global + (local_num_cells - 1)

          ! Set number of vertices per cell
          call set_vert_per_cell(4_ccs_int, mesh)
          call set_vert_nb_per_cell(4_ccs_int, mesh)

          call get_vert_per_cell(mesh, vert_per_cell)
          call get_vert_nb_per_cell(mesh, vert_nb_per_cell)

          ! Allocate mesh topolgy arrays
          allocate (mesh%topo%global_indices(local_num_cells))
          allocate (mesh%topo%num_nb(local_num_cells))
          allocate (mesh%topo%num_vert_nb(local_num_cells))
          allocate (mesh%topo%nb_indices(max_faces, local_num_cells))
          allocate (mesh%topo%vert_nb_indices(vert_nb_per_cell, local_num_cells))
          allocate (mesh%topo%face_indices(max_faces, local_num_cells))

          ! Initialise mesh arrays
          mesh%topo%num_nb(:) = max_faces ! All cells have 4 neighbours (possibly ghost/boundary cells)
          mesh%topo%num_vert_nb(:) = vert_nb_per_cell ! All cells have 4 vertex neighbours (possibly ghost/boundary cells)

          ! Initialise neighbour indices
          mesh%topo%nb_indices(:, :) = 0_ccs_int

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
            nb_direction = 0_ccs_int
            set_vert_nb = .false.

            ! Construct left (1) face/neighbour
            nb_direction(1) = left
            face_counter = left
            call add_neighbour(i, face_counter, index_counter, nb_direction, cps, cps, cps, set_vert_nb, mesh)

            ! Construct right (2) face/neighbour
            nb_direction(1) = right
            face_counter = right
            call add_neighbour(i, face_counter, index_counter, nb_direction, cps, cps, cps, set_vert_nb, mesh)

            ! Construct bottom (3) face/neighbour
            nb_direction(1) = bottom
            face_counter = bottom
            call add_neighbour(i, face_counter, index_counter, nb_direction, cps, cps, cps, set_vert_nb, mesh)

            ! Construct top (4) face/neighbour
            nb_direction(1) = top
            face_counter = top
            call add_neighbour(i, face_counter, index_counter, nb_direction, cps, cps, cps, set_vert_nb, mesh)

            ! Now construct vertex neighbours
            set_vert_nb = .true.
            nb_direction = [top, left]
            vertex_counter = front_top_left
            call add_neighbour(i, vertex_counter, index_counter, nb_direction, cps, cps, cps, set_vert_nb, mesh)

            ! Construct top right neighbour
            nb_direction = [top, right]
            vertex_counter = front_top_right
            call add_neighbour(i, vertex_counter, index_counter, nb_direction, cps, cps, cps, set_vert_nb, mesh)

            ! Construct bottom left neighbour
            nb_direction = [bottom, left]
            vertex_counter = front_bottom_left
            call add_neighbour(i, vertex_counter, index_counter, nb_direction, cps, cps, cps, set_vert_nb, mesh)

            ! Construct bottom right neighbour
            nb_direction = [bottom, right]
            vertex_counter = front_bottom_right
            call add_neighbour(i, vertex_counter, index_counter, nb_direction, cps, cps, cps, set_vert_nb, mesh)

            index_counter = index_counter + 1_ccs_int
          end do
        end associate

        call set_total_num_cells(size(mesh%topo%global_indices), mesh)
        call get_total_num_cells(mesh, total_num_cells)
        call set_halo_num_cells(total_num_cells - local_num_cells, mesh)

        ! Create shared memory global arrays
        call create_shared_array(shared_env, global_num_faces, mesh%topo%bnd_rid, mesh%topo%bnd_rid_window)

        length(1) = max_faces
        length(2) = nglobal
        call create_shared_array(shared_env, length(:), mesh%topo%global_face_indices, mesh%topo%global_face_indices_window)

        ! Initialise shared memory global arrays
        if (is_root(shared_env)) then
          mesh%topo%global_face_indices(:, :) = 0_ccs_int
          mesh%topo%bnd_rid(:) = 0_ccs_int
        end if

        ! Construct face_cell1 and face_cell2 following:
        !  - face_cell1 < face_cell2
        !  - and if face is a boundary, then: face_cell1 = current_cell, face_cell2 = 0
        face_index_counter = 1_ccs_int
        do i = 1, nglobal

          ii = i - 1_ccs_int

          ! Construct left (1) face/neighbour
          face_counter = left
          if (modulo(ii, cps) == 0_ccs_int) then
            global_index_nb = -left

            call create_face_locator(mesh, i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
            face_index_counter = face_index_counter + 1_ccs_int
          else
            ! If internal left face, nothing to be done, the face will be linked as a right face from its neighbour
            !global_index_nb = i - 1_ccs_int
            !mesh%topo%face_cell1(face_index_counter) = global_index_nb
            !mesh%topo%face_cell2(face_index_counter) = i
          end if

          ! Construct right (2) face/neighbour
          face_counter = right
          if (modulo(ii, cps) == (cps - 1_ccs_int)) then
            global_index_nb = -right

            call create_face_locator(mesh, i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
          else
            global_index_nb = i + 1_ccs_int

            call create_face_locator(mesh, i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            call create_face_locator(mesh, global_index_nb, left, loc_f)
            call set_global_index(face_index_counter, loc_f)
          end if
          face_index_counter = face_index_counter + 1_ccs_int

          ! Construct bottom (3) face/neighbour
          face_counter = bottom
          if (modulo(ii / cps, cps) == 0_ccs_int) then
            global_index_nb = -bottom

            call create_face_locator(mesh, i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
            face_index_counter = face_index_counter + 1_ccs_int
          else
            ! If internal bottom face, nothing to be done, the face will be linked as a top face from its neighbour
            !global_index_nb = i - nx
            !mesh%topo%face_cell1(face_index_counter) = global_index_nb
            !mesh%topo%face_cell2(face_index_counter) = i
          end if

          ! Construct top (4) face/neighbour
          face_counter = top
          if (modulo(ii / cps, cps) == (cps - 1_ccs_int)) then
            global_index_nb = -top

            call create_face_locator(mesh, i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
          else
            global_index_nb = i + cps

            call create_face_locator(mesh, i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            call create_face_locator(mesh, global_index_nb, bottom, loc_f)
            call set_global_index(face_index_counter, loc_f)
          end if
          face_index_counter = face_index_counter + 1_ccs_int

        end do

        call set_num_faces(count_mesh_faces(mesh), mesh)

        call set_cell_face_indices(mesh)

        length(1) = mesh%topo%vert_per_cell
        length(2) = mesh%topo%global_num_cells
       call create_shared_array(shared_env, length, mesh%topo%global_vertex_indices, mesh%topo%global_vertex_indices_window)

        ! Global vertex numbering
        if (is_root(shared_env)) then
          do i = 1, mesh%topo%global_num_cells
            ii = i
            associate (global_vert_index => mesh%topo%global_vertex_indices(:, i))

              global_vert_index(front_bottom_left) = ii + (ii - 1) / cps
              global_vert_index(front_bottom_right) = global_vert_index(front_bottom_left) + 1
              global_vert_index(front_top_left) = global_vert_index(front_bottom_left) + (cps + 1)
              global_vert_index(front_top_right) = global_vert_index(front_top_left) + 1
            end associate
          end do
        end if
        call sync(shared_env)

      class default
        call error_abort("Unknown parallel environment type.")

      end select

    class default
      call error_abort("Unknown parallel environment type.")

    end select

  end subroutine build_square_topology

  subroutine build_square_topology_connectivity(shared_env, &
                                                cps, &
                                                global_num_faces, max_faces, &
                                                face_cell1, face_cell1_window, face_cell2, face_cell2_window)

    class(parallel_environment), intent(in) :: shared_env !< The shared parallel environment
    integer(ccs_int), intent(in) :: cps
    integer(ccs_int), intent(out) :: global_num_faces
    integer(ccs_int), intent(out) :: max_faces
    integer(ccs_int), dimension(:), pointer, intent(out) :: face_cell1
    integer, intent(out) :: face_cell1_window
    integer(ccs_int), dimension(:), pointer, intent(out) :: face_cell2
    integer, intent(out) :: face_cell2_window

    integer(ccs_int) :: i
    integer(ccs_int) :: ii
    integer(ccs_int) :: nglobal
    integer(ccs_int) :: face_index_counter
    integer(ccs_int) :: global_index_nb

    nglobal = cps**2
    global_num_faces = 2 * cps * (cps + 1)
    max_faces = 4 ! Constant for square meshes

    call create_shared_array(shared_env, global_num_faces, face_cell1, face_cell1_window)
    call create_shared_array(shared_env, global_num_faces, face_cell2, face_cell2_window)

    if (is_root(shared_env)) then
      face_cell1(:) = 0_ccs_int
      face_cell2(:) = 0_ccs_int
    end if

    face_index_counter = 1
    do i = 1, nglobal
      ii = i - 1

      ! Left face
      if (modulo(ii, cps) == 0_ccs_int) then
        face_cell1(face_index_counter) = i
        face_cell2(face_index_counter) = 0
        face_index_counter = face_index_counter + 1_ccs_int
      end if

      ! Right face
      if (modulo(ii, cps) == (cps - 1_ccs_int)) then
        face_cell1(face_index_counter) = i
        face_cell2(face_index_counter) = 0
      else
        global_index_nb = i + 1_ccs_int
        face_cell1(face_index_counter) = i
        face_cell2(face_index_counter) = global_index_nb
      end if
      face_index_counter = face_index_counter + 1_ccs_int

      ! Bottom face
      if (modulo(ii / cps, cps) == 0_ccs_int) then
        face_cell1(face_index_counter) = i
        face_cell2(face_index_counter) = 0
        face_index_counter = face_index_counter + 1_ccs_int
      end if

      ! Top face
      if (modulo(ii / cps, cps) == (cps - 1_ccs_int)) then
        face_cell1(face_index_counter) = i
        face_cell2(face_index_counter) = 0
      else
        global_index_nb = i + cps
        face_cell1(face_index_counter) = i
        face_cell2(face_index_counter) = global_index_nb
      end if
      face_index_counter = face_index_counter + 1_ccs_int

    end do

  end subroutine build_square_topology_connectivity

  subroutine build_square_geometry(par_env, cps, side_length, mesh)

    class(parallel_environment), intent(in) :: par_env !< The parallel environment to construct the mesh.
    integer(ccs_int), intent(in) :: cps                !< Number of cells per side of the mesh.
    real(ccs_real), intent(in) :: side_length          !< The length of each side.

    type(ccs_mesh), intent(inout) :: mesh                             !< The resulting mesh.

    integer(ccs_int) :: i               ! Loop counter
    integer(ccs_int) :: ii              ! Zero-indexed loop counter (simplifies some operations)
    integer(ccs_int) :: face_counter    ! Cell-local face counter
    integer(ccs_int) :: vertex_counter  ! Cell-local vertex counter
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: total_num_cells
    integer(ccs_int) :: max_faces
    integer(ccs_int) :: vert_per_cell

    logical :: is_boundary

    integer(ccs_int) :: index_nb        ! The local index of a neighbour cell

    real(ccs_real), dimension(2) :: x_p ! Cell centre array
    type(cell_locator) :: loc_p         ! Cell locator object

    real(ccs_real), dimension(3) :: x_nb_3 ! Cell centre array of neighbour cell
    real(ccs_real), dimension(2) :: x_nb   ! Cell centre array of neighbour cell
    type(neighbour_locator) :: loc_nb    ! the neighbour locator object.

    real(ccs_real), dimension(2) :: x_f    ! Face centre array
    real(ccs_real), dimension(2) :: normal ! Face normal array
    type(face_locator) :: loc_f            ! Face locator object

    real(ccs_real), dimension(2) :: x_v ! Vertex centre array
    type(vert_locator) :: loc_v         ! Vertex locator object

    select type (par_env)
    type is (parallel_environment_mpi)

      call get_local_num_cells(mesh, local_num_cells)
      call get_vert_per_cell(mesh, vert_per_cell)

      call get_total_num_cells(mesh, total_num_cells)
      call get_max_faces(mesh, max_faces)
      allocate (mesh%geo%x_p(ndim, total_num_cells))
      allocate (mesh%geo%x_f(ndim, max_faces, local_num_cells)) !< @note Currently hardcoded as a 2D mesh. @endnote
      allocate (mesh%geo%volumes(total_num_cells))
      allocate (mesh%geo%face_areas(max_faces, local_num_cells))
      allocate (mesh%geo%face_normals(ndim, max_faces, local_num_cells)) ! Currently hardcoded as a 2D mesh.
      allocate (mesh%geo%vert_coords(ndim, vert_per_cell, local_num_cells))

      mesh%geo%h = side_length / real(cps, ccs_real)
      mesh%geo%volumes(:) = mesh%geo%h**2 !< @note Mesh is square and 2D @endnote
      mesh%geo%face_normals(:, :, :) = 0.0_ccs_real
      mesh%geo%x_p(:, :) = 0.0_ccs_real
      mesh%geo%x_f(:, :, :) = 0.0_ccs_real
      mesh%geo%face_areas(:, :) = mesh%geo%h  ! Mesh is square and 2D
      mesh%geo%vert_coords(:, :, :) = 0.0_ccs_real

      ! Set cell centre
      associate (h => mesh%geo%h)
        do i = 1_ccs_int, total_num_cells
          call create_cell_locator(mesh, i, loc_p)
          call get_natural_index(loc_p, ii)

          x_p(1) = (modulo(ii - 1, cps) + 0.5_ccs_real) * h
          x_p(2) = ((ii - 1) / cps + 0.5_ccs_real) * h

          call set_centre(loc_p, x_p)
        end do

        do i = 1_ccs_int, local_num_cells
          call create_cell_locator(mesh, i, loc_p)
          call get_centre(loc_p, x_p)

          do face_counter = 1_ccs_int, max_faces

            call create_neighbour_locator(loc_p, face_counter, loc_nb)
            call get_boundary_status(loc_nb, is_boundary)
            call create_face_locator(mesh, i, face_counter, loc_f)

            if (.not. is_boundary) then
              ! faces are midway between cell centre and nb cell centre
              call get_centre(loc_nb, x_nb_3)
              x_nb(:) = x_nb_3(1:2) ! XXX: hacky fix for issue with resolving get_neighbour_centre in 2D

              x_f(:) = 0.5_ccs_real * (x_p(:) + x_nb(:))
              normal(:) = (x_nb(:) - x_p(:)) / h
              call set_centre(loc_f, x_f)
              call set_normal(loc_f, normal)

            else
              ! for boundary faces we use their 'ID' to get their location

              call get_local_index(loc_nb, index_nb)
              x_f(1) = x_p(1)
              x_f(2) = x_p(2)
              normal(1) = 0.0_ccs_real
              normal(2) = 0.0_ccs_real

              if (index_nb .eq. -left) then
                x_f(1) = x_p(1) - 0.5_ccs_real * h
                normal(1) = -1.0_ccs_real
              end if

              if (index_nb .eq. -right) then
                x_f(1) = x_p(1) + 0.5_ccs_real * h
                normal(1) = +1.0_ccs_real
              end if

              if (index_nb .eq. -bottom) then
                x_f(2) = x_p(2) - 0.5_ccs_real * h
                normal(2) = -1.0_ccs_real
              end if

              if (index_nb .eq. -top) then
                x_f(2) = x_p(2) + 0.5_ccs_real * h
                normal(2) = +1.0_ccs_real
              end if

              call set_centre(loc_f, x_f)
              call set_normal(loc_f, normal)

            end if
          end do
        end do

        do i = 1_ccs_int, local_num_cells
          call create_cell_locator(mesh, i, loc_p)
          call get_centre(loc_p, x_p)

          vertex_counter = front_bottom_left
          call create_vert_locator(mesh, i, vertex_counter, loc_v)
          x_v(1) = x_p(1) - 0.5_ccs_real * h
          x_v(2) = x_p(2) - 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = front_bottom_right
          call create_vert_locator(mesh, i, vertex_counter, loc_v)
          x_v(1) = x_p(1) + 0.5_ccs_real * h
          x_v(2) = x_p(2) - 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = front_top_left
          call create_vert_locator(mesh, i, vertex_counter, loc_v)
          x_v(1) = x_p(1) - 0.5_ccs_real * h
          x_v(2) = x_p(2) + 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = front_top_right
          call create_vert_locator(mesh, i, vertex_counter, loc_v)
          x_v(1) = x_p(1) + 0.5_ccs_real * h
          x_v(2) = x_p(2) + 0.5_ccs_real * h
          call set_centre(loc_v, x_v)
        end do
      end associate

      call compute_face_interpolation(mesh)

    class default
      call error_abort("Unknown parallel environment type.")

    end select

  end subroutine build_square_geometry

  !v Utility constructor to build a 3D mesh with hex cells.
  !
  !  Builds a Cartesian grid of nx*ny*nz cells.
  function build_mesh(par_env, shared_env, nx, ny, nz, side_length) result(mesh)

    use partitioning, only: compute_partitioner_input

    class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared memory environment
    integer(ccs_int), intent(in) :: nx                 !< Number of cells in the x direction.
    integer(ccs_int), intent(in) :: ny                 !< Number of cells in the y direction.
    integer(ccs_int), intent(in) :: nz                 !< Number of cells in the z direction.
    real(ccs_real), intent(in) :: side_length          !< The length of the side.

    type(ccs_mesh) :: mesh                             !< The resulting mesh.

    character(:), allocatable :: error_message
    call set_mesh_generated(.true., mesh)

    if (.not. (nx .eq. ny .and. ny .eq. nz)) then !< @note Must be a cube (for now) @endnote
      error_message = "Only supporting cubes for now - nx, ny and nz must be the same!"
      call error_abort(error_message)
    end if

    if (nx * ny * nz < par_env%num_procs) then
      error_message = "ERROR: Global number of cells < number of ranks. &
                      &Increase the mesh size or reduce the number of MPI ranks."
      call error_abort(error_message)
    end if

    call build_topology(par_env, shared_env, nx, ny, nz, mesh)

    call compute_partitioner_input(par_env, shared_env, mesh)

    call mesh_partition_reorder(par_env, shared_env, mesh)

    call build_geometry(par_env, nx, ny, nz, side_length, mesh)

    call cleanup_topo(shared_env, mesh)

  end function build_mesh

  !v Utility constructor to build a 3D mesh with hex cells.
  !
  !  Builds a Cartesian grid of nx*ny*nz cells.
  subroutine build_topology(par_env, shared_env, nx, ny, nz, mesh)

    class(parallel_environment), intent(in) :: par_env    !< The parallel environment to construct the mesh.
    class(parallel_environment), intent(in) :: shared_env !< The shared memory environment
    integer(ccs_int), intent(in) :: nx                 !< Number of cells in the x direction.
    integer(ccs_int), intent(in) :: ny                 !< Number of cells in the y direction.
    integer(ccs_int), intent(in) :: nz                 !< Number of cells in the z direction.

    type(ccs_mesh), intent(inout) :: mesh                             !< The resulting mesh.

    integer(ccs_int) :: start_global    ! The (global) starting index of a partition
    integer(ccs_int) :: end_global      ! The (global) last index of a partition
    integer(ccs_int) :: i               ! Loop counter
    integer(ccs_int) :: ii              ! Zero-indexed loop counter (simplifies some operations)
    integer(ccs_int) :: index_counter   ! Local index counter
    integer(ccs_int) :: face_counter    ! Cell-local face counter
    integer(ccs_int) :: face_index_counter    ! global face counter
    integer(ccs_int) :: vertex_counter  ! Cell-local vertex counter
    integer(ccs_int) :: a, b, c, d, e       ! Temporary variables

    integer(ccs_int) :: global_index_nb ! The global index of a neighbour cell
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: total_num_cells
    integer(ccs_int) :: global_num_faces
    integer(ccs_int) :: max_faces
    integer(ccs_int) :: vert_per_cell
    integer(ccs_int) :: vert_nb_per_cell
    integer(ccs_int) :: global_num_vertices
    integer(ccs_int) :: nglobal
    logical :: set_vert_nb              ! Flag for setting vertex neighbour
    integer(ccs_int), dimension(3) :: nb_direction  ! Array indicating direction of neighbour

    type(face_locator) :: loc_f

    integer(ccs_int), dimension(2) :: length

    nglobal = nx * ny * nz ! The global cell count
    call build_topology_connectivity(shared_env, &
                                     nx, ny, nz, &
                                     global_num_faces, max_faces, &
                                     mesh%topo%face_cell1, mesh%topo%face_cell1_window, &
                                     mesh%topo%face_cell2, mesh%topo%face_cell2_window)
    call set_naive_distribution(par_env, nglobal, mesh%topo%graph_conn)

    !< XXX: <It should be possible to enter the partitioner here>

    select type (par_env)
    type is (parallel_environment_mpi)

      select type (shared_env)
      type is (parallel_environment_mpi)

        ! Set the global mesh parameters
        call set_global_num_cells(nglobal, mesh)
        call set_global_num_vertices((nx + 1) * (ny + 1) * (nz + 1), mesh)
        call set_global_num_faces(global_num_faces, mesh)
        call set_max_faces(max_faces, mesh)

        ! Just to make sure we are working with the same numbers as the mesh object.
        call get_global_num_cells(mesh, nglobal)
        call get_global_num_vertices(mesh, global_num_vertices)
        call get_global_num_faces(mesh, global_num_faces)
        call get_max_faces(mesh, max_faces)

        ! Determine ownership range
        start_global = global_start(nglobal, par_env%proc_id, par_env%num_procs)
        local_num_cells = local_count(nglobal, par_env%proc_id, par_env%num_procs)
        call set_local_num_cells(local_num_cells, mesh)
        call get_local_num_cells(mesh, local_num_cells) ! Ensure using correct value

        ! Abort the execution if any rank has 0 local cells
        if (local_num_cells <= 0) then
          call error_abort("ERROR: Zero local cells found.")
        end if

        call set_total_num_cells(local_num_cells, mesh) ! Setting initial value
        end_global = start_global + (local_num_cells - 1)

        ! Set max number of faces (constant, 6)
        call set_max_faces(6_ccs_int, mesh)

        ! Set number of vertices per cell (constant, 8)
        call set_vert_per_cell(8, mesh)

        ! Set number of neighbours via vertex per cell
        call set_vert_nb_per_cell(20_ccs_int, mesh)

        call get_max_faces(mesh, max_faces)
        call get_vert_per_cell(mesh, vert_per_cell)
        call get_vert_nb_per_cell(mesh, vert_nb_per_cell)

        ! Allocate mesh arrays
        allocate (mesh%topo%global_indices(local_num_cells))
        allocate (mesh%topo%num_nb(local_num_cells))
        allocate (mesh%topo%num_vert_nb(local_num_cells))
        allocate (mesh%topo%nb_indices(max_faces, local_num_cells))
        allocate (mesh%topo%vert_nb_indices(vert_nb_per_cell, local_num_cells))
        allocate (mesh%topo%face_indices(max_faces, local_num_cells))

        ! Initialise mesh arrays
        mesh%topo%num_nb(:) = max_faces ! All cells have 6 neighbours (possibly ghost/boundary cells)
        mesh%topo%num_vert_nb(:) = vert_nb_per_cell

        ! Initalise neighbour indices
        mesh%topo%nb_indices(:, :) = 0_ccs_int
        mesh%topo%vert_nb_indices(:, :) = 0_ccs_int

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
          set_vert_nb = .false.
          nb_direction(:) = 0_ccs_int

          ! Construct left (1) face/neighbour
          nb_direction(1) = left
          face_counter = left
          call add_neighbour(i, face_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          ! Construct right (2) face/neighbour
          nb_direction(1) = right
          face_counter = right
          call add_neighbour(i, face_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          ! Construct bottom (3) face/neighbour
          nb_direction(1) = bottom
          face_counter = bottom
          call add_neighbour(i, face_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          ! Construct top (4) face/neighbour
          nb_direction(1) = top
          face_counter = top
          call add_neighbour(i, face_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          ! Construct back (5) face/neighbour
          nb_direction(1) = back
          face_counter = back
          call add_neighbour(i, face_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          ! Construct front (6) face/neighbour
          nb_direction(1) = front
          face_counter = front
          call add_neighbour(i, face_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          ! Now construct neighbours connected via vertex or edge.
          ! There are 8 front neighbours, 4 middle neighbours and 8 back neighbours
          set_vert_nb = .true.
          nb_direction = [front, top, left]
          vertex_counter = front_top_left
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          nb_direction = [front, top, 0_ccs_int]
          vertex_counter = front_top
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          nb_direction = [front, top, right]
          vertex_counter = front_top_right
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          nb_direction = [front, right, 0_ccs_int]
          vertex_counter = front_right
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          nb_direction = [front, bottom, right]
          vertex_counter = front_bottom_right
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          nb_direction = [front, bottom, 0_ccs_int]
          vertex_counter = front_bottom
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          nb_direction = [front, bottom, left]
          vertex_counter = front_bottom_left
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          nb_direction = [front, left, 0_ccs_int]
          vertex_counter = front_left
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          ! Now do the middle layer
          nb_direction = [top, left, 0_ccs_int]
          vertex_counter = middle_top_left
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          nb_direction = [top, right, 0_ccs_int]
          vertex_counter = middle_top_right
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          nb_direction = [bottom, right, 0_ccs_int]
          vertex_counter = middle_bottom_right
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          nb_direction = [bottom, left, 0_ccs_int]
          vertex_counter = middle_bottom_left
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          ! And finally the back layer, again start at top left
          nb_direction = [back, top, left]
          vertex_counter = back_top_left
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          nb_direction = [back, top, 0_ccs_int]
          vertex_counter = back_top
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          nb_direction = [back, top, right]
          vertex_counter = back_top_right
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          nb_direction = [back, right, 0_ccs_int]
          vertex_counter = back_right
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          nb_direction = [back, bottom, right]
          vertex_counter = back_bottom_right
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          nb_direction = [back, bottom, 0_ccs_int]
          vertex_counter = back_bottom
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          nb_direction = [back, bottom, left]
          vertex_counter = back_bottom_left
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          nb_direction = [back, left, 0_ccs_int]
          vertex_counter = back_left
          call add_neighbour(i, vertex_counter, index_counter, nb_direction, nx, ny, nz, set_vert_nb, mesh)

          index_counter = index_counter + 1_ccs_int

        end do

        ! print*,"Neighbour indices: ",mesh%neighbour_indices

        call set_total_num_cells(size(mesh%topo%global_indices), mesh)
        call get_total_num_cells(mesh, total_num_cells)
        call set_halo_num_cells(total_num_cells - local_num_cells, mesh)

        ! Create shared memory global arrays
        call create_shared_array(shared_env, global_num_faces, mesh%topo%bnd_rid, mesh%topo%bnd_rid_window)

        length(1) = max_faces
        length(2) = nglobal
        call create_shared_array(shared_env, length(:), mesh%topo%global_face_indices, mesh%topo%global_face_indices_window)

        ! Initialise shared memory global arrays
        if (is_root(shared_env)) then
          mesh%topo%global_face_indices(:, :) = 0_ccs_int
          mesh%topo%bnd_rid(:) = 0_ccs_int
        end if

        ! Construct face_cell1 and face_cell2 following:
        !  - face_cell1 < face_cell2
        !  - and if face is a boundary, then: face_cell1 = current_cell, face_cell2 = 0
        face_index_counter = 1_ccs_int

        do i = 1, nglobal

          ii = i - 1_ccs_int

          ! Construct left (1) face/neighbour
          face_counter = left
          if (modulo(ii, nx) == 0_ccs_int) then
            global_index_nb = -left

            call create_face_locator(mesh, i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
            face_index_counter = face_index_counter + 1_ccs_int
          else
            ! If internal left face, nothing to be done, the face will be linked as a right face from its neighbour
            !global_index_nb = i - 1_ccs_int
            !mesh%topo%face_cell1(face_index_counter) = global_index_nb
            !mesh%topo%face_cell2(face_index_counter) = i
          end if

          ! Construct right (2) face/neighbour
          face_counter = right
          if (modulo(ii, nx) == (nx - 1_ccs_int)) then
            global_index_nb = -right

            call create_face_locator(mesh, i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
          else
            global_index_nb = i + 1_ccs_int

            call create_face_locator(mesh, i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            call create_face_locator(mesh, global_index_nb, left, loc_f)
            call set_global_index(face_index_counter, loc_f)
          end if
          face_index_counter = face_index_counter + 1_ccs_int

          ! Construct bottom (3) face/neighbour
          face_counter = bottom
          if (modulo(ii / nx, ny) == 0_ccs_int) then
            global_index_nb = -bottom

            call create_face_locator(mesh, i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
            face_index_counter = face_index_counter + 1_ccs_int
          else
            ! If internal bottom face, nothing to be done, the face will be linked as a top face from its neighbour
            !global_index_nb = i - nx
            !mesh%topo%face_cell1(face_index_counter) = global_index_nb
            !mesh%topo%face_cell2(face_index_counter) = i
          end if

          ! Construct top (4) face/neighbour
          face_counter = top
          if (modulo(ii / nx, ny) == (ny - 1_ccs_int)) then
            global_index_nb = -top

            call create_face_locator(mesh, i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
          else
            global_index_nb = i + nx

            call create_face_locator(mesh, i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            call create_face_locator(mesh, global_index_nb, bottom, loc_f)
            call set_global_index(face_index_counter, loc_f)
          end if
          face_index_counter = face_index_counter + 1_ccs_int

          ! Construct back (5) face/neighbour
          face_counter = back
          if ((ii / (nx * ny)) == 0_ccs_int) then
            global_index_nb = -back

            call create_face_locator(mesh, i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
            face_index_counter = face_index_counter + 1_ccs_int
          else
            ! If internal back face, nothing to be done, the face will be linked as a front face from its neighbour
            !global_index_nb = i - nx * ny
            !mesh%topo%face_cell1(face_index_counter) = global_index_nb
            !mesh%topo%face_cell2(face_index_counter) = i
          end if

          ! Construct front (6) face/neighbour
          face_counter = front
          if ((ii / (nx * ny)) == nz - 1_ccs_int) then
            global_index_nb = -front

            call create_face_locator(mesh, i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
          else
            global_index_nb = i + nx * ny

            call create_face_locator(mesh, i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            call create_face_locator(mesh, global_index_nb, back, loc_f)
            call set_global_index(face_index_counter, loc_f)
          end if
          face_index_counter = face_index_counter + 1_ccs_int

        end do

        call set_num_faces(count_mesh_faces(mesh), mesh)

        call set_cell_face_indices(mesh)

        length(1) = mesh%topo%vert_per_cell
        length(2) = mesh%topo%global_num_cells
       call create_shared_array(shared_env, length, mesh%topo%global_vertex_indices, mesh%topo%global_vertex_indices_window)

        ! Global vertex numbering
        if (is_root(shared_env)) then
          do i = 1, mesh%topo%global_num_cells
            associate (global_vert_index => mesh%topo%global_vertex_indices(:, i))
              ii = i
              a = modulo(ii - 1, nx * ny) + 1
              b = (a - 1) / nx
              c = ((ii - 1) / (nx * ny)) * (nx + 1) * (ny + 1)
              d = (a + nx - 1) / nx
              e = (nx + 1) * (ny + 1)

              global_vert_index(front_bottom_left) = a + b + c
              global_vert_index(front_bottom_right) = a + b + c + 1
              global_vert_index(front_top_left) = a + c + d + nx
              global_vert_index(front_top_right) = a + c + d + nx + 1
              global_vert_index(back_bottom_left) = a + b + c + e
              global_vert_index(back_bottom_right) = a + b + c + e + 1
              global_vert_index(back_top_left) = a + c + d + e + nx
              global_vert_index(back_top_right) = a + c + d + e + nx + 1
            end associate
          end do
        end if
        call sync(shared_env)

      class default
        call error_abort("Unknown parallel environment type.")

      end select

    class default
      call error_abort("Unknown parallel environment type.")

    end select

  end subroutine build_topology

  subroutine build_topology_connectivity(shared_env, &
                                         nx, ny, nz, &
                                         global_num_faces, max_faces, &
                                         face_cell1, face_cell1_window, face_cell2, face_cell2_window)

    class(parallel_environment), intent(in) :: shared_env !< The shared parallel environment
    integer(ccs_int), intent(in) :: nx, ny, nz
    integer(ccs_int), intent(out) :: global_num_faces
    integer(ccs_int), intent(out) :: max_faces
    integer(ccs_int), dimension(:), pointer, intent(out) :: face_cell1
    integer, intent(out) :: face_cell1_window
    integer(ccs_int), dimension(:), pointer, intent(out) :: face_cell2
    integer, intent(out) :: face_cell2_window

    integer(ccs_int) :: i
    integer(ccs_int) :: ii
    integer(ccs_int) :: nglobal
    integer(ccs_int) :: face_index_counter
    integer(ccs_int) :: global_index_nb

    nglobal = nx * ny * nz
    global_num_faces = (nx + 1) * ny * nz + nx * (ny + 1) * nz + nx * ny * (nz + 1)
    max_faces = 6 ! Constant for hex meshes

    call create_shared_array(shared_env, global_num_faces, face_cell1, face_cell1_window)
    call create_shared_array(shared_env, global_num_faces, face_cell2, face_cell2_window)

    if (is_root(shared_env)) then
      face_cell1(:) = 0_ccs_int
      face_cell2(:) = 0_ccs_int
    end if

    face_index_counter = 1
    do i = 1, nglobal
      ii = i - 1

      ! Left face
      if (modulo(ii, nx) == 0_ccs_int) then
        face_cell1(face_index_counter) = i
        face_cell2(face_index_counter) = 0
        face_index_counter = face_index_counter + 1_ccs_int
      end if

      ! Right face
      if (modulo(ii, nx) == (nx - 1_ccs_int)) then
        face_cell1(face_index_counter) = i
        face_cell2(face_index_counter) = 0
      else
        global_index_nb = i + 1_ccs_int
        face_cell1(face_index_counter) = i
        face_cell2(face_index_counter) = global_index_nb
      end if
      face_index_counter = face_index_counter + 1_ccs_int

      ! Bottom face
      if (modulo(ii / nx, ny) == 0_ccs_int) then
        face_cell1(face_index_counter) = i
        face_cell2(face_index_counter) = 0
        face_index_counter = face_index_counter + 1_ccs_int
      end if

      ! Top face
      if (modulo(ii / nx, ny) == (ny - 1_ccs_int)) then
        face_cell1(face_index_counter) = i
        face_cell2(face_index_counter) = 0
      else
        global_index_nb = i + nx
        face_cell1(face_index_counter) = i
        face_cell2(face_index_counter) = global_index_nb
      end if
      face_index_counter = face_index_counter + 1_ccs_int

      ! Back face
      if ((ii / (nx * ny)) == 0_ccs_int) then
        face_cell1(face_index_counter) = i
        face_cell2(face_index_counter) = 0
        face_index_counter = face_index_counter + 1_ccs_int
      end if

      ! Front face
      if ((ii / (nx * ny)) == nz - 1_ccs_int) then
        global_index_nb = -front
        face_cell1(face_index_counter) = i
        face_cell2(face_index_counter) = 0
      else
        global_index_nb = i + nx * ny
        face_cell1(face_index_counter) = i
        face_cell2(face_index_counter) = global_index_nb
      end if
      face_index_counter = face_index_counter + 1_ccs_int

    end do

  end subroutine build_topology_connectivity

  !v Utility constructor to build a 3D mesh with hex cells.
  !
  !  Builds a Cartesian grid of nx*ny*nz cells.
  subroutine build_geometry(par_env, nx, ny, nz, side_length, mesh)

    class(parallel_environment), intent(in) :: par_env !< The parallel environment to construct the mesh.
    integer(ccs_int), intent(in) :: nx                 !< Number of cells in the x direction.
    integer(ccs_int), intent(in) :: ny                 !< Number of cells in the y direction.
    integer(ccs_int), intent(in) :: nz                 !< Number of cells in the z direction.
    real(ccs_real), intent(in) :: side_length          !< The length of the side.

    type(ccs_mesh), intent(inout) :: mesh                             !< The resulting mesh.

    integer(ccs_int) :: i               ! Loop counter
    integer(ccs_int) :: ii              ! Zero-indexed loop counter (simplifies some operations)
    integer(ccs_int) :: face_counter    ! Cell-local face counter
    integer(ccs_int) :: vertex_counter  ! Cell-local vertex counter

    logical :: is_boundary

    integer(ccs_int) :: index_nb        ! The local index of a neighbour cell
    integer(ccs_int) :: local_num_cells ! The local cell count
    integer(ccs_int) :: total_num_cells ! The total cell count
    integer(ccs_int) :: max_faces       ! The maximum number of faces per cell
    integer(ccs_int) :: vert_per_cell

    real(ccs_real), dimension(3) :: x_p ! Cell centre array
    type(cell_locator) :: loc_p         ! Cell locator object

    real(ccs_real), dimension(3) :: x_nb ! Cell centre array of neighbour cell
    type(neighbour_locator) :: loc_nb    ! the neighbour locator object.

    real(ccs_real), dimension(3) :: x_f    ! Face centre array
    real(ccs_real), dimension(3) :: normal ! Face normal array
    type(face_locator) :: loc_f            ! Face locator object

    real(ccs_real), dimension(3) :: x_v ! Vertex centre array
    type(vert_locator) :: loc_v         ! Vertex locator object

    associate (foo => nz) ! Silence unused dummy argument
    end associate

    call get_local_num_cells(mesh, local_num_cells) ! Ensure using correct value
    call get_vert_per_cell(mesh, vert_per_cell)

    select type (par_env)
    type is (parallel_environment_mpi)

      call get_total_num_cells(mesh, total_num_cells)
      call get_max_faces(mesh, max_faces)
      allocate (mesh%geo%x_p(ndim, total_num_cells))
      allocate (mesh%geo%x_f(ndim, max_faces, local_num_cells))
      allocate (mesh%geo%volumes(total_num_cells))
      allocate (mesh%geo%face_areas(max_faces, local_num_cells))
      allocate (mesh%geo%face_normals(ndim, max_faces, local_num_cells))
      allocate (mesh%geo%vert_coords(ndim, vert_per_cell, local_num_cells))

      mesh%geo%h = side_length / real(nx, ccs_real) !< @note Assumes cube @endnote
      mesh%geo%volumes(:) = mesh%geo%h**3 !< @note Mesh is cube @endnote
      mesh%geo%face_normals(:, :, :) = 0.0_ccs_real
      mesh%geo%x_p(:, :) = 0.0_ccs_real
      mesh%geo%x_f(:, :, :) = 0.0_ccs_real
      mesh%geo%face_areas(:, :) = mesh%geo%h**2
      mesh%geo%vert_coords(:, :, :) = 0.0_ccs_real

      ! Set cell centre
      associate (h => mesh%geo%h)
        do i = 1_ccs_int, total_num_cells
          call create_cell_locator(mesh, i, loc_p)
          call get_natural_index(loc_p, ii)

          x_p(1) = (modulo(ii - 1, nx) + 0.5_ccs_real) * h
          x_p(2) = (modulo((ii - 1) / nx, ny) + 0.5_ccs_real) * h
          x_p(3) = (((ii - 1) / (nx * ny)) + 0.5_ccs_real) * h

          call set_centre(loc_p, x_p)
        end do

        do i = 1_ccs_int, local_num_cells
          call create_cell_locator(mesh, i, loc_p)
          call get_centre(loc_p, x_p)

          do face_counter = 1_ccs_int, max_faces

            call create_neighbour_locator(loc_p, face_counter, loc_nb)
            call get_boundary_status(loc_nb, is_boundary)
            call create_face_locator(mesh, i, face_counter, loc_f)

            if (.not. is_boundary) then
              ! faces are midway between cell centre and nb cell centre
              call get_centre(loc_nb, x_nb)

              x_f(:) = 0.5_ccs_real * (x_p(:) + x_nb(:))
              normal(:) = (x_nb(:) - x_p(:)) / h
              call set_centre(loc_f, x_f)
              call set_normal(loc_f, normal)

            else
              ! for boundary faces we use their 'ID' to get their location

              call get_local_index(loc_nb, index_nb)
              x_f(1) = x_p(1)
              x_f(2) = x_p(2)
              x_f(3) = x_p(3)
              normal(1) = 0.0_ccs_real
              normal(2) = 0.0_ccs_real
              normal(3) = 0.0_ccs_real

              if (index_nb .eq. -left) then
                x_f(1) = x_p(1) - 0.5_ccs_real * h
                normal(1) = -1.0_ccs_real
              end if

              if (index_nb .eq. -right) then
                x_f(1) = x_p(1) + 0.5_ccs_real * h
                normal(1) = +1.0_ccs_real
              end if

              if (index_nb .eq. -bottom) then
                x_f(2) = x_p(2) - 0.5_ccs_real * h
                normal(2) = -1.0_ccs_real
              end if

              if (index_nb .eq. -top) then
                x_f(2) = x_p(2) + 0.5_ccs_real * h
                normal(2) = +1.0_ccs_real
              end if

              if (index_nb .eq. -back) then
                x_f(3) = x_p(3) - 0.5_ccs_real * h
                normal(3) = -1.0_ccs_real
              end if

              if (index_nb .eq. -front) then
                x_f(3) = x_p(3) + 0.5_ccs_real * h
                normal(3) = +1.0_ccs_real
              end if

              call set_centre(loc_f, x_f)
              call set_normal(loc_f, normal)

            end if

          end do
        end do

        do i = 1_ccs_int, local_num_cells
          call create_cell_locator(mesh, i, loc_p)
          call get_centre(loc_p, x_p)

          vertex_counter = front_bottom_left
          call create_vert_locator(mesh, i, vertex_counter, loc_v)
          x_v(1) = x_p(1) - 0.5_ccs_real * h
          x_v(2) = x_p(2) - 0.5_ccs_real * h
          x_v(3) = x_p(3) + 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = front_bottom_right
          call create_vert_locator(mesh, i, vertex_counter, loc_v)
          x_v(1) = x_p(1) + 0.5_ccs_real * h
          x_v(2) = x_p(2) - 0.5_ccs_real * h
          x_v(3) = x_p(3) + 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = front_top_left
          call create_vert_locator(mesh, i, vertex_counter, loc_v)
          x_v(1) = x_p(1) - 0.5_ccs_real * h
          x_v(2) = x_p(2) + 0.5_ccs_real * h
          x_v(3) = x_p(3) + 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = front_top_right
          call create_vert_locator(mesh, i, vertex_counter, loc_v)
          x_v(1) = x_p(1) + 0.5_ccs_real * h
          x_v(2) = x_p(2) + 0.5_ccs_real * h
          x_v(3) = x_p(3) + 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = back_bottom_left
          call create_vert_locator(mesh, i, vertex_counter, loc_v)
          x_v(1) = x_p(1) - 0.5_ccs_real * h
          x_v(2) = x_p(2) - 0.5_ccs_real * h
          x_v(3) = x_p(3) - 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = back_bottom_right
          call create_vert_locator(mesh, i, vertex_counter, loc_v)
          x_v(1) = x_p(1) + 0.5_ccs_real * h
          x_v(2) = x_p(2) - 0.5_ccs_real * h
          x_v(3) = x_p(3) - 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = back_top_left
          call create_vert_locator(mesh, i, vertex_counter, loc_v)
          x_v(1) = x_p(1) - 0.5_ccs_real * h
          x_v(2) = x_p(2) + 0.5_ccs_real * h
          x_v(3) = x_p(3) - 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = back_top_right
          call create_vert_locator(mesh, i, vertex_counter, loc_v)
          x_v(1) = x_p(1) + 0.5_ccs_real * h
          x_v(2) = x_p(2) + 0.5_ccs_real * h
          x_v(3) = x_p(3) - 0.5_ccs_real * h
          call set_centre(loc_v, x_v)
        end do
      end associate

      call compute_face_interpolation(mesh)

    class default
      call error_abort("Unknown parallel environment type.")

    end select

  end subroutine build_geometry

  !> Helper subroutine to appropriately set local and global neighbour indices
  subroutine add_neighbour(index_p, nb_counter, index_counter, direction, nx, ny, nz, vertex_flag, mesh)
    integer(ccs_int), intent(in) :: index_p                   !< Global index of cell whose neighbours we're assembling
    integer(ccs_int), intent(in) :: nb_counter                !< the cell-relative index neighbour index
    integer(ccs_int), intent(in) :: index_counter             !< local index of cell whose neighbours we're assembling
    integer(ccs_int), dimension(:), intent(in) :: direction   !< Array containing the direction of the neighbour
    !< relative to the cell
    integer(ccs_int), intent(in) :: nx                        !< Mesh size in x direction
    integer(ccs_int), intent(in) :: ny                        !< Mesh size in y direction
    integer(ccs_int), intent(in) :: nz                        !< Mesh size in z direction
    logical, intent(in) :: vertex_flag                        !< Flag to indicate whether this is a vertex neighbour
    type(ccs_mesh), intent(inout) :: mesh                     !< The mesh

    integer(ccs_int) :: i, ii
    integer(ccs_int) :: index_nb
    integer(ccs_int) :: global_index_nb
    integer(ccs_int) :: index_increment

    index_nb = 0_ccs_int

    ii = index_p - 1    ! We're indexing cells starting from 1, so adjust for modulo computations.
    index_increment = 0_ccs_int

    do i = 1, size(direction)
      select case (direction(i))
      case (left)
        if (modulo(ii, nx) == 0_ccs_int) then
          index_nb = -left
          global_index_nb = -left
        end if
        index_increment = index_increment - 1_ccs_int
      case (right)
        if (modulo(ii, nx) == (nx - 1_ccs_int)) then
          index_nb = -right
          global_index_nb = -right
        end if
        index_increment = index_increment + 1_ccs_int
      case (bottom)
        if (modulo(ii / nx, ny) == 0_ccs_int) then
          index_nb = -bottom
          global_index_nb = -bottom
        end if
        index_increment = index_increment - nx
      case (top)
        if (modulo(ii / nx, ny) == (ny - 1_ccs_int)) then
          index_nb = -top
          global_index_nb = -top
        end if
        index_increment = index_increment + nx
      case (back)
        if ((ii / (nx * ny)) == 0_ccs_int) then
          index_nb = -back
          global_index_nb = -back
        end if
        index_increment = index_increment - nx * ny
      case (front)
        if ((ii / (nx * ny)) == nz - 1_ccs_int) then
          index_nb = -front
          global_index_nb = -front
        end if
        index_increment = index_increment + nx * ny
      case default
        if (direction(i) /= 0) then
          call error_abort("Unexpected value for neighbour direction. Index " // str(i) // " direction " // str(direction(i)))
        end if
      end select
    end do

    if (index_nb == 0) then
      index_nb = index_counter + index_increment
      global_index_nb = index_p + index_increment
    end if

    call build_local_mesh_add_neighbour(index_counter, nb_counter, index_nb, global_index_nb, mesh, vertex_flag)
  end subroutine add_neighbour

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
  subroutine build_local_mesh_add_neighbour(index_p, index_p_nb, index_nb, global_index_nb, mesh, vertex_nb_flag)

    integer(ccs_int), intent(in) :: index_p !< the index of the cell whose neighbours we are assembling
    integer(ccs_int), intent(in) :: index_p_nb !< the cell-relative neighbour index
    integer(ccs_int), intent(in) :: index_nb !< the local index of the neighbour cell
    integer(ccs_int), intent(in) :: global_index_nb !< the global index of the neighbour cell
    type(ccs_mesh), intent(inout) :: mesh !< the mesh we are assembling neighbours on
    logical, intent(in) :: vertex_nb_flag !< flag indicating whether the neighbour being added is a vertex neighbour

    integer(ccs_int) :: ng  ! The current number of cells (total = local + halos)
    logical :: found        ! Indicates whether a halo cell was already present
    integer(ccs_int) :: i   ! Cell iteration counter

    integer(ccs_int) :: local_num_cells ! The number of local cells
    integer(ccs_int) :: global_num_cells ! The number of global cells

    type(cell_locator) :: loc_h
    integer(ccs_int) :: global_index_h

    integer(ccs_int) :: total_num_cells

    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb

    call get_local_num_cells(mesh, local_num_cells)
    call get_global_num_cells(mesh, global_num_cells)

    call create_cell_locator(mesh, index_p, loc_p)
    if (.not. vertex_nb_flag) then
      call create_neighbour_locator(loc_p, index_p_nb, loc_nb)
    end if
    if ((index_nb >= 1_ccs_int) .and. (index_nb <= local_num_cells)) then
      ! Neighbour is local
      if (vertex_nb_flag) then
        mesh%topo%vert_nb_indices(index_p_nb, index_p) = index_nb
      else
        call set_local_index(index_nb, loc_nb)
      end if
    else if (global_index_nb < 0_ccs_int) then
      ! Boundary "neighbour" - local index should also be -ve
      if (.not. (index_nb < 0_ccs_int)) then
        call error_abort("ERROR: boundary neighbours should have -ve indices.")
      end if
      if (vertex_nb_flag) then
        mesh%topo%vert_nb_indices(index_p_nb, index_p) = index_nb
      else
        call set_local_index(index_nb, loc_nb)
      end if
    else
      ! Neighbour is in a halo

      ! First check if neighbour is already present in halo
      ng = size(mesh%topo%global_indices)
      found = .false.
      do i = local_num_cells + 1, ng
        call create_cell_locator(mesh, i, loc_h)
        call get_global_index(loc_h, global_index_h)
        if (global_index_h == global_index_nb) then
          found = .true.
          if (vertex_nb_flag) then
            mesh%topo%vert_nb_indices(index_p_nb, index_p) = i
          else
            call set_local_index(i, loc_nb)
          end if
          exit
        end if
      end do

      ! If neighbour was not present append to global index list (the end of the global index list
      ! becoming its local index).
      ! XXX: Note this currently copies into an n+1 temporary, reallocates and then copies back to
      !      the (extended) original array.
      if (.not. found) then
        if ((ng + 1) > global_num_cells) then
          call error_abort("ERROR: Trying to create halo that exceeds global mesh size.")
        end if

        call append_to_arr(global_index_nb, mesh%topo%global_indices)
        ng = size(mesh%topo%global_indices)
        if (vertex_nb_flag) then
          mesh%topo%vert_nb_indices(index_p_nb, index_p) = ng
        else
          call set_local_index(ng, loc_nb)
        end if

        ! Increment total cell count
        call get_total_num_cells(mesh, total_num_cells)
        call set_total_num_cells(total_num_cells + 1, mesh)

        call get_total_num_cells(mesh, total_num_cells)
        if (total_num_cells /= size(mesh%topo%global_indices)) then
          print *, total_num_cells, size(mesh%topo%global_indices)
          call error_abort("ERROR: Local total cell count and size of global indices not in agreement")
        end if
      end if
    end if

  end subroutine build_local_mesh_add_neighbour

  !v @note Docs needed.
  subroutine append_to_arr(i, arr)

    integer(ccs_int), intent(in) :: i
    integer(ccs_int), dimension(:), allocatable, intent(inout) :: arr ! XXX: Allocatable here be
    ! dragons. If this were intent(out) it
    ! would be deallocated on entry
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
    integer(ccs_int) :: n_faces_internal ! Internal face count
    integer(ccs_int) :: nfaces_bnd       ! Boundary face count
    integer(ccs_int) :: nfaces_interface ! Process interface face count
    logical :: is_boundary
    logical :: is_local
    integer(ccs_int) :: local_num_cells

    ! Initialise
    n_faces_internal = 0
    nfaces_bnd = 0
    nfaces_interface = 0

    ! Loop over cells
    call get_local_num_cells(mesh, local_num_cells)
    do index_p = 1, local_num_cells
      call create_cell_locator(mesh, index_p, loc_p)
      call get_global_index(loc_p, global_index_p)
      call count_neighbours(loc_p, nnb)

      do j = 1, nnb
        call create_neighbour_locator(loc_p, j, loc_nb)
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
    integer(ccs_int) :: face_counter      ! Face index counter
    logical :: is_boundary
    integer(ccs_int) :: local_num_cells

    face_counter = 0

    ! Loop over cells
    call get_local_num_cells(mesh, local_num_cells)
    do index_p = 1, local_num_cells
      call create_cell_locator(mesh, index_p, loc_p)
      call count_neighbours(loc_p, nnb)

      do j = 1, nnb
        call create_neighbour_locator(loc_p, j, loc_nb)
        call get_local_index(loc_nb, index_nb)
        call get_boundary_status(loc_nb, is_boundary)

        if (.not. is_boundary) then
          ! Cell with lowest local index assigns an index to the face
          if (index_p < index_nb) then
            face_counter = face_counter + 1
            call set_face_index(index_p, j, face_counter, mesh)
          else
            ! Find corresponding face in neighbour cell
            ! (To be improved, this seems inefficient)
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

    call create_cell_locator(mesh, index_nb, loc_nb)
    call count_neighbours(loc_nb, nnb_nb)
    do k = 1, nnb_nb
      call create_neighbour_locator(loc_nb, k, loc_nb_nb)
      call get_local_index(loc_nb_nb, index_nb_nb)
      if (index_nb_nb == index_p) then
        call create_face_locator(mesh, index_nb, k, loc_f)
        call get_local_index(loc_f, index_f)
        exit ! Exit the loop, as found shared face
      else if (k == nnb_nb) then
        call error_abort("ERROR: Failed to find face in owning cell.")
      end if
    end do
  end function get_neighbour_face_index

  !v From cell centre and face centre, computes face interpolation factors
  ! face interpolations are bounded between 0 and 1, and are relative to
  ! the cell with the lowest cell index.
  subroutine compute_face_interpolation(mesh)
    type(ccs_mesh), intent(inout) :: mesh

    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb
    type(face_locator) :: loc_f

    real(ccs_real) :: interpol_factor
    integer(ccs_int) :: index_p, index_nb
    integer(ccs_int) :: j
    integer(ccs_int) :: nnb
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: num_faces
    logical :: is_boundary

    real(ccs_real), dimension(ndim) :: x_p ! cell centre array
    real(ccs_real), dimension(ndim) :: x_nb ! neighbour cell centre array
    real(ccs_real), dimension(ndim) :: x_f ! face centre array
    real(ccs_real), dimension(ndim) :: v_p_nb ! vector going through local and neighbour cell centres

    if (allocated(mesh%geo%face_interpol)) then
      deallocate (mesh%geo%face_interpol)
    end if

    call get_num_faces(mesh, num_faces)
    allocate (mesh%geo%face_interpol(num_faces))

    ! Safe guard to make sure we go through all faces
    mesh%geo%face_interpol(:) = -1.0_ccs_real

    call get_local_num_cells(mesh, local_num_cells)

    do index_p = 1, local_num_cells
      call create_cell_locator(mesh, index_p, loc_p)
      call count_neighbours(loc_p, nnb)

      do j = 1, nnb
        call create_neighbour_locator(loc_p, j, loc_nb)
        call get_boundary_status(loc_nb, is_boundary)
        call create_face_locator(mesh, index_p, j, loc_f)

        if (.not. is_boundary) then
          call get_local_index(loc_nb, index_nb)

          call get_centre(loc_f, x_f)
          call get_centre(loc_nb, x_nb)
          call get_centre(loc_p, x_p)

          ! v_p_nb.V2 / |v_p_nb|**2
          v_p_nb = x_nb - x_p
          interpol_factor = dot_product(v_p_nb, x_f - x_p) / dot_product(v_p_nb, v_p_nb)

          ! inverse interpol factor as it is relative to x_p
          ! the closer x_f is to x_p, the higher the interpol_factor
          interpol_factor = 1.0_ccs_real - interpol_factor

          call set_face_interpolation(interpol_factor, loc_f)

        else
          ! Boundary faces values are not meaningful and shouldn't be read
          call set_face_interpolation(0.0_ccs_real, loc_f)
        end if
      end do  ! End loop over current cell's neighbours
    end do    ! End loop over local cells

    if (minval(mesh%geo%face_interpol) < 0.0_ccs_real .or. &
        maxval(mesh%geo%face_interpol) > 1.0_ccs_real) then
      call error_abort("Face interpolation out of bound.")
    end if

  end subroutine

  ! Populate mesh%topo%graph_conn%global_partition with a split of cells in stride using global_start and local_count
  subroutine partition_stride(par_env, shared_env, roots_env, mesh)
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: roots_env !< The parallel environment
    type(ccs_mesh), intent(inout) :: mesh                             !< The resulting mesh.

    integer(ccs_int) :: iproc, start, end
    integer(ccs_int) :: global_num_cells

    ! roots_env kept as argument for consistency with partition_kway
    associate (foo => roots_env)
    end associate

    call get_global_num_cells(mesh, global_num_cells)

    if (is_root(shared_env)) then
      do iproc = 0, par_env%num_procs - 1
        start = global_start(global_num_cells, iproc, par_env%num_procs)
        end = start + local_count(global_num_cells, iproc, par_env%num_procs) - 1
        mesh%topo%graph_conn%global_partition(start:end) = iproc
      end do
    end if

  end subroutine partition_stride

  !v Compute the global start index of a process if each receives an equal share of the items,
  !  ordered by process ID. Any remainder items are distributed evenly between the remainder lower
  !  process IDs.
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

  !v Compute the local share of a process if each receives an equal share of the items, ordered by
  !  process ID. Any remainder items are distributed evenly between the remainder lower process IDs.
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

  ! Print mesh geometry object
  subroutine print_geo(par_env, mesh)
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(ccs_mesh), intent(in) :: mesh       !< the mesh
    integer(ccs_int) :: i          ! loop counters
    integer(ccs_int) :: nb_elem = 10

    print *, par_env%proc_id, "############################# Print Geometry ########################################"

    print *, par_env%proc_id, "h                  : ", mesh%geo%h
    print *, par_env%proc_id, "scalefactor        : ", mesh%geo%scalefactor
    print *, ""

    if (allocated(mesh%geo%volumes)) then
      print *, par_env%proc_id, "volumes     : ", mesh%geo%volumes(1:nb_elem)
    else
      print *, par_env%proc_id, "volumes     : UNALLOCATED"
    end if

    if (allocated(mesh%geo%face_interpol)) then
      print *, par_env%proc_id, "face_interpol          : ", mesh%geo%face_interpol(1:nb_elem)
    else
      print *, par_env%proc_id, "face_interpol          : UNALLOCATED"
    end if

    print *, ""
    if (allocated(mesh%geo%face_areas)) then
      do i = 1, nb_elem
        print *, par_env%proc_id, "face_areas(1:" // str(nb_elem / 2) // ", " // str(i) // ")", &
          mesh%geo%face_areas(1:nb_elem / 2, i)
      end do
    else
      print *, par_env%proc_id, "face_areas             : UNALLOCATED"
    end if

    print *, ""
    if (allocated(mesh%geo%x_p)) then
      do i = 1, nb_elem
        print *, par_env%proc_id, "x_p(:)", mesh%geo%x_p(:, i)
      end do
    else
      print *, par_env%proc_id, "x_p                    : UNALLOCATED"
    end if

    print *, ""
    if (allocated(mesh%geo%x_f)) then
      do i = 1, nb_elem
        print *, par_env%proc_id, "x_f(2, 1:" // str(nb_elem / 2) // ", " // str(i) // ")", &
          mesh%geo%x_f(2, 1:nb_elem / 2, i)
      end do
    else
      print *, par_env%proc_id, "x_f                    : UNALLOCATED"
    end if

    print *, ""
    if (allocated(mesh%geo%face_normals)) then
      do i = 1, nb_elem
        print *, par_env%proc_id, "face_normals(2, 1:" // str(nb_elem / 2) // ", " // str(i) // ")", &
          mesh%geo%face_normals(2, 1:nb_elem / 2, i)
      end do
    else
      print *, par_env%proc_id, "face_normals          : UNALLOCATED"
    end if

    print *, ""
    if (allocated(mesh%geo%vert_coords)) then
      do i = 1, nb_elem
        print *, par_env%proc_id, "vert_coords(2, 1:" // str(nb_elem / 2) // ", " // str(i) // ")", &
          mesh%geo%vert_coords(2, 1:nb_elem / 2, i)
      end do
    else
      print *, par_env%proc_id, "vert_coords           : UNALLOCATED"
    end if

    print *, par_env%proc_id, "############################# End Print Geometry ########################################"

  end subroutine print_geo

  ! Print mesh topology object
  subroutine print_topo(par_env, mesh)

    use meshing, only: get_halo_num_cells

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(ccs_mesh), intent(in) :: mesh       !< the mesh
    integer(ccs_int) :: i          ! loop counters
    integer(ccs_int) :: nb_elem = 10

    integer(ccs_int) :: local_num_cells, global_num_cells, halo_num_cells, total_num_cells
    integer(ccs_int) :: global_num_faces, num_faces, max_faces
    integer(ccs_int) :: vert_per_cell, global_num_vertices

    print *, par_env%proc_id, "############################# Print Topology ########################################"

    call get_local_num_cells(mesh, local_num_cells)
    call get_global_num_cells(mesh, global_num_cells)
    call get_halo_num_cells(mesh, halo_num_cells)
    call get_total_num_cells(mesh, total_num_cells)
    call get_global_num_faces(mesh, global_num_faces)
    call get_num_faces(mesh, num_faces)
    call get_max_faces(mesh, max_faces)
    call get_vert_per_cell(mesh, vert_per_cell)
    call get_global_num_vertices(mesh, global_num_vertices)

    print *, par_env%proc_id, "global_num_cells    : ", global_num_cells
    print *, par_env%proc_id, "local_num_cells     : ", local_num_cells
    print *, par_env%proc_id, "halo_num_cells      : ", halo_num_cells
    print *, par_env%proc_id, "total_num_cells     : ", total_num_cells
    print *, par_env%proc_id, "global_num_vertices : ", global_num_vertices
    print *, par_env%proc_id, "vert_per_cell       : ", vert_per_cell
    print *, par_env%proc_id, "global_num_faces    : ", global_num_faces
    print *, par_env%proc_id, "num_faces           : ", num_faces
    print *, par_env%proc_id, "max_faces           : ", max_faces
    print *, ""

    if (allocated(mesh%topo%global_indices)) then
      print *, par_env%proc_id, "global_indices     : ", mesh%topo%global_indices(1:nb_elem)
    else
      print *, par_env%proc_id, "global_indices     : UNALLOCATED"
    end if

    if (allocated(mesh%topo%num_nb)) then
      print *, par_env%proc_id, "num_nb             : ", mesh%topo%num_nb(1:nb_elem)
    else
      print *, par_env%proc_id, "num_nb             : UNALLOCATED"
    end if

    if (associated(mesh%topo%face_cell1)) then
      print *, par_env%proc_id, "face_cell1        : ", mesh%topo%face_cell1(1:nb_elem)
    else
      print *, par_env%proc_id, "face_cell1        : UNALLOCATED"
    end if

    if (associated(mesh%topo%face_cell2)) then
      print *, par_env%proc_id, "face_cell2        : ", mesh%topo%face_cell2(1:nb_elem)
    else
      print *, par_env%proc_id, "face_cell2        : UNALLOCATED"
    end if

    if (associated(mesh%topo%bnd_rid)) then
      print *, par_env%proc_id, "bnd_rid           : ", mesh%topo%bnd_rid(1:nb_elem)
    else
      print *, par_env%proc_id, "bnd_rid           : UNALLOCATED"
    end if

    if (allocated(mesh%topo%graph_conn%vwgt)) then
      print *, par_env%proc_id, "vwgt              : ", mesh%topo%graph_conn%vwgt(1:nb_elem)
    else
      print *, par_env%proc_id, "vwgt              : UNALLOCATED"
    end if

    if (allocated(mesh%topo%graph_conn%adjwgt)) then
      print *, par_env%proc_id, "adjwgt            : ", mesh%topo%graph_conn%adjwgt(1:nb_elem)
    else
      print *, par_env%proc_id, "adjwgt            : UNALLOCATED"
    end if

    if (allocated(mesh%topo%graph_conn%local_partition)) then
      print *, par_env%proc_id, "local_partition   : ", mesh%topo%graph_conn%local_partition(1:nb_elem)
    else
      print *, par_env%proc_id, "local_partition   : UNALLOCATED"
    end if

    if (associated(mesh%topo%graph_conn%global_partition)) then
      print *, par_env%proc_id, "global_partition  : ", mesh%topo%graph_conn%global_partition(1:nb_elem)
    else
      print *, par_env%proc_id, "global_partition  : UNALLOCATED"
    end if

    print *, ""
    if (associated(mesh%topo%global_face_indices)) then
      do i = 1, nb_elem
        print *, par_env%proc_id, "global_face_indices(1:"  //  str(nb_elem / 2)  //  ", "  //  str(i)  //  ")", mesh%topo%global_face_indices(1:nb_elem / 2, i)
      end do
    else
      print *, par_env%proc_id, "global_face_indices   : UNALLOCATED"
    end if

    print *, ""
    if (associated(mesh%topo%global_vertex_indices)) then
      do i = 1, nb_elem
        print *, par_env%proc_id, "global_vertex_indices(1:"  //  str(nb_elem / 2)  //  ", "  //  str(i)  //  ")", mesh%topo%global_vertex_indices(1:nb_elem / 2, i)
      end do
    else
      print *, par_env%proc_id, "global_vertex_indices : UNALLOCATED"
    end if

    print *, ""
    if (allocated(mesh%topo%face_indices)) then
      do i = 1, nb_elem
    print *, par_env%proc_id, "face_indices(1:" // str(nb_elem / 2) // ", " // str(i) // ")", mesh%topo%face_indices(1:nb_elem / 2, i)
      end do
    else
      print *, par_env%proc_id, "face_indices          : UNALLOCATED"
    end if

    print *, ""
    if (allocated(mesh%topo%nb_indices)) then
      do i = 1, nb_elem
        print *, par_env%proc_id, "nb_indices(1:" // str(nb_elem / 2) // ", " // str(i) // ")", mesh%topo%nb_indices(1:nb_elem / 2, i)
      end do
    else
      print *, par_env%proc_id, "nb_indices            : UNALLOCATED"
    end if

    print *, par_env%proc_id, "############################# End Print Topology ########################################"

  end subroutine print_topo

  subroutine cleanup_topo(shared_env, mesh)

    type(ccs_mesh), target, intent(inout) :: mesh   !< The mesh
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The parallel environment

    if (associated(mesh%topo%global_face_indices)) then
      call destroy_shared_array(shared_env, mesh%topo%global_face_indices, mesh%topo%global_face_indices_window)
      call dprint("mesh%topo%global_face_indices deallocated.")
    end if

    if (associated(mesh%topo%global_vertex_indices)) then
      call destroy_shared_array(shared_env, mesh%topo%global_vertex_indices, mesh%topo%global_vertex_indices_window)
      call dprint("mesh%topo%global_vertex_indices deallocated.")
    end if

    if (associated(mesh%topo%face_cell1)) then
      call destroy_shared_array(shared_env, mesh%topo%face_cell1, mesh%topo%face_cell1_window)
      call dprint("mesh%topo%face_cell1 deallocated.")
    end if

    if (associated(mesh%topo%face_cell2)) then
      call destroy_shared_array(shared_env, mesh%topo%face_cell2, mesh%topo%face_cell2_window)
      call dprint("mesh%topo%face_cell2 deallocated.")
    end if

    if (associated(mesh%topo%bnd_rid)) then
      call destroy_shared_array(shared_env, mesh%topo%bnd_rid, mesh%topo%bnd_rid_window)
      call dprint("mesh%topo%bnd_rid deallocated.")
    end if

    if (associated(mesh%topo%graph_conn%global_partition)) then
      call destroy_shared_array(shared_env, mesh%topo%graph_conn%global_partition, mesh%topo%graph_conn%global_partition_window)
      call dprint("mesh%topo%graph_conn%global_partition deallocated.")
    end if
  end subroutine cleanup_topo

  subroutine mesh_partition_reorder(par_env, shared_env, mesh)

    use partitioning, only: partition_kway, &
                            compute_connectivity, &
                            compute_partitioner_input, &
                            cleanup_partitioner_data
    use case_config, only: compute_bwidth

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The parallel environment
    class(parallel_environment), allocatable, target :: roots_env !< The parallel environment
    type(ccs_mesh), intent(inout) :: mesh                                   !< The mesh

    integer(ccs_int) :: bw_max
    real(ccs_real) :: bw_avg

    call create_shared_roots_comm(par_env, shared_env, roots_env)

    if (par_env%num_procs > 1) then
      call partition_kway(par_env, shared_env, roots_env, mesh)
    else
      call partition_stride(par_env, shared_env, roots_env, mesh)
    end if

    call compute_connectivity(par_env, shared_env, roots_env, mesh)

    if(compute_bwidth .eqv. .true.) then
      call compute_bandwidth(mesh, bw_max, bw_avg)
      call print_bandwidth(par_env, bw_max, bw_avg)
    end if

    call reorder_cells(par_env, shared_env, mesh)
    call cleanup_partitioner_data(shared_env, mesh)

    if(compute_bwidth .eqv. .true.) then
      call compute_bandwidth(mesh, bw_max, bw_avg)
      call print_bandwidth(par_env, bw_max, bw_avg)
    end if

  end subroutine

  subroutine print_bandwidth(par_env, bw_max, bw_avg)
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    integer(ccs_int), intent(in) :: bw_max
    real(ccs_real), intent(in) :: bw_avg

    integer(ccs_int) :: ierr
    real(ccs_real) :: sum_bw_avg

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


  ! Naively distribute cells equally across all processes
  subroutine set_naive_distribution(par_env, num_cells, graph_conn)

    class(parallel_environment), intent(in) :: par_env
    integer(ccs_int), intent(in) :: num_cells
    type(graph_connectivity), intent(inout) :: graph_conn

    integer(ccs_int) :: i, j, k

    ! Create and populate the vtxdist array based on the total number of cells
    ! and the total number of ranks in the parallel environment
    allocate (graph_conn%vtxdist(par_env%num_procs + 1)) ! vtxdist array is of size num_procs + 1 on all ranks

    graph_conn%vtxdist(1) = 1                                        ! First element is 1
    graph_conn%vtxdist(par_env%num_procs + 1) = num_cells + 1 ! Last element is total number of cells + 1

    ! Divide the total number of cells by the world size to compute the chunk sizes
    k = int(real(num_cells) / par_env%num_procs)
    j = 1

    do i = 1, par_env%num_procs
      graph_conn%vtxdist(i) = j
      j = j + k
    end do
  end subroutine set_naive_distribution

end module mesh_utils
