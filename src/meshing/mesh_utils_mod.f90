!v Module file mesh_utils_mod.f90
!
!  Implements utilities for loading/generating meshes.
!
!  @dont_fail_linter
!  This file triggers many false alarms in the linter.

module mesh_utils
#include "ccs_macros.inc"

  use mpi

  use profiler
  use core, only: ccs_options
  use constants, only: ndim, geoext, adiosconfig
  use utils, only: exit_print, str, debug_print
  use kinds, only: ccs_int, ccs_long, ccs_real, ccs_err
  use types, only: ccs_mesh, topology, &
                   io_environment, io_process, &
                   face_locator, cell_locator, neighbour_locator, vert_locator, &
                   graph_connectivity
  use io, only: read_scalar, read_array, &
                configure_io, open_file, close_file, &
                initialise_io, cleanup_io
  use parallel, only: create_shared_array, is_root, is_valid, create_shared_roots_comm, &
                      destroy_shared_array, sync
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi
  use meshing, only: get_global_index, get_natural_index, get_local_index, count_neighbours, &
                     create_cell_locator, create_neighbour_locator, create_face_locator, create_vert_locator, &
                     set_face_index, get_boundary_status, get_local_status, &
                     get_centre, set_centre, &
                     set_area, set_normal, get_face_normal, get_face_area, &
                     set_total_num_cells, get_total_num_cells, &
                     get_local_num_cells, set_local_num_cells, set_face_interpolation, &
                     get_global_num_cells, set_global_num_cells, &
                     set_halo_num_cells, &
                     get_global_num_faces, set_global_num_faces, &
                     get_num_faces, set_num_faces, &
                     get_max_faces, set_max_faces, &
                     get_vert_per_cell, set_vert_per_cell, &
                     get_global_num_vertices, set_global_num_vertices, &
                     set_face_interpolation, &
                     set_local_index, &
                     set_global_index, &
                     set_mesh_generated, set_mesh_object, nullify_mesh_object, &
                     set_topo_object, nullify_topo_object
  use bc_constants
  use reordering, only: reorder_cells, print_bandwidth
  use logging, only: log_unit_out

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
  public :: destroy_mesh_storage
  public :: global_start
  public :: local_count
  public :: count_mesh_faces
  public :: set_cell_face_indices
  public :: compute_face_interpolation
  public :: partition_stride
  public :: print_topo
  public :: print_geo
  public :: build_adjacency_matrix
  public :: test_mesh_internal_neighbours

contains

  !v Read mesh from file
  subroutine read_mesh(par_env, shared_env, run_options, mesh)

    use partitioning, only: compute_partitioner_input

    use profiler

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The parallel environment
    type(ccs_options), intent(in) :: run_options
    type(ccs_mesh), intent(inout) :: mesh                                   !< The mesh

    ! Local variables
    character(len=:), allocatable :: geo_file    ! Geo file name
    character(len=:), allocatable :: adios2_file ! ADIOS2 config file name

    class(parallel_environment), allocatable, target :: reader_env !< The reader parallel environment
    class(io_environment), allocatable :: io_env
    class(io_process), allocatable :: geo_reader

    character(len=:), allocatable :: case_name

    call set_mesh_object(mesh)
    call set_mesh_generated(.false.)
    call nullify_mesh_object()

    case_name = run_options%paths%case_name
    geo_file = case_name // "_mesh" // geoext
    adios2_file = case_name // adiosconfig
    if (is_root(par_env)) then
      write (log_unit_out, *) "Mesh file:", geo_file
    end if

    call create_shared_roots_comm(par_env, shared_env, reader_env)

    call initialise_io(reader_env, adios2_file, io_env)
    call configure_io(io_env, "geo_reader", geo_reader)

    call open_file(geo_file, "read", geo_reader)

    call profiler_begin_region("Read mesh topology")
    call read_topology(par_env, shared_env, reader_env, geo_reader, mesh)
    call profiler_end_region("Read mesh topology")

    call profiler_begin_region("Compute partitioner input")
    call compute_partitioner_input(par_env, shared_env, mesh)
    call profiler_end_region("Compute partitioner input")

    call mesh_partition_reorder(par_env, shared_env, run_options, mesh)

    call set_offsets(shared_env, mesh)

    call profiler_begin_region("Read mesh geometry")
    call read_geometry(shared_env, reader_env, geo_reader, mesh)
    call profiler_end_region("Read mesh geometry")

    ! Close the file and ADIOS2 engine
    call close_file(geo_reader)

    ! Finalise the ADIOS2 IO environment
    call cleanup_io(io_env)

    ! TODO: cleanup reader parallel environment

    call cleanup_topo(shared_env, mesh)

    mesh%bnd_names = run_options%mesh%bnd_names
    call check_mesh_bnd_names(par_env, mesh)

    ! Compute the mesh surface integral if the mesh was read from an input file
    call profiler_begin_region("Compute mesh surface integral")
    call report_mesh_surface_integral(par_env, mesh)
    call profiler_end_region("Compute mesh surface integral")

  end subroutine read_mesh

  !> Helper subroutine to check boundary ID/name compatibility.
  subroutine check_mesh_bnd_names(par_env, mesh)

    class(parallel_environment), intent(in) :: par_env
    type(ccs_mesh), intent(in) :: mesh

    integer :: i
    integer :: bc_cnt

    logical :: id_names_valid

    integer :: ierr

    if (is_root(par_env)) then
      write (log_unit_out, *) "=========================="
      write (log_unit_out, *) "Boundary ID map"
      do i = 1, size(mesh%bnd_names)
        write (log_unit_out, *) i, trim(mesh%bnd_names(i))
      end do
    end if

    ! The most negative value of neighbour indices (i.e. maximum boundary ID) should equal the
    ! boundary count. Note that in parallel any given process may not have this many boundaries so
    ! we use a logical OR to check at least one process has a valid boundary ID count.
    ! We need to also check that no process exceeds the boundary ID count.
    bc_cnt = -minval(mesh%topo%nb_indices)

    ! Check range
    id_names_valid = (bc_cnt == size(mesh%bnd_names))
    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_Allreduce(MPI_IN_PLACE, id_names_valid, 1, MPI_LOGICAL, MPI_LOR, par_env%comm, ierr)
    class default
      call error_abort("Unsupported parallel environment")
    end select
    if (.not. id_names_valid) then
      call error_abort("Maximum boundary ID doesn't match supplied boundary name count")
    end if

    ! Check no boundary IDs exceed the range
    id_names_valid = (bc_cnt <= size(mesh%bnd_names))
    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_Allreduce(MPI_IN_PLACE, id_names_valid, 1, MPI_LOGICAL, MPI_LOR, par_env%comm, ierr)
    class default
      call error_abort("Unsupported parallel environment")
    end select
    if (.not. id_names_valid) then
      call error_abort("Maximum boundary ID doesn't match supplied boundary name count")
    end if
    call sync(par_env)

    if (is_root(par_env)) then
      write (log_unit_out, *) "Boundary name list / ID compatibility: PASS"
      write (log_unit_out, *) "=========================="
    end if

  end subroutine check_mesh_bnd_names

  !v Read the topology data from an input (HDF5) file
  ! This subroutine assumes the following names are used in the file:
  ! "ncel" - the total number of cells
  ! "nfac" - the total number of faces
  ! "maxfaces" - the maximum number of faces per cell
  ! "/face/cell1" and "/face/cell2" - the arrays the face edge data
  !
  ! This high-level interface zeroes the topology object contained by the mesh before calling the
  ! lower-level routine to read the topology object.
  subroutine read_topology(par_env, shared_env, reader_env, geo_reader, mesh)
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: reader_env !< The reader parallel environment
    class(io_process), intent(in) :: geo_reader                                         !< The IO process for reading the file
    type(ccs_mesh), intent(inout) :: mesh                                   !< The mesh that will be read

    call set_mesh_object(mesh)

    ! Zero scalar topology values to have known initial state
    call set_global_num_cells(0_ccs_int)
    call set_global_num_faces(0_ccs_int)
    call set_max_faces(0_ccs_int)
    call set_local_num_cells(0_ccs_int)
    call set_halo_num_cells(0_ccs_int)
    call set_total_num_cells(0_ccs_int)
    call set_num_faces(0_ccs_int)
    call set_vert_per_cell(0_ccs_int)
    call set_global_num_vertices(0_ccs_int)

    call nullify_mesh_object()

    call read_topology_topo(par_env, shared_env, reader_env, geo_reader, mesh%topo)

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
    class(io_process), intent(in) :: geo_reader                                         !< The IO process for reading the file
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

    call set_topo_object(topo)
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
    call get_global_num_cells(global_num_cells)

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
    call get_max_faces(max_faces)
    if (max_faces == 6) then ! if cell are hexes
      call set_vert_per_cell(8) ! 8 vertices per cell
    else if (max_faces == 4) then ! if cell are tetrahedral
      call set_vert_per_cell(4) ! 4 vertices per cell
    else
      call error_abort("Currently only supporting hex or tet cells.")
    end if

    call get_global_num_faces(global_num_faces)

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
    call get_vert_per_cell(vert_per_cell)
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
      call set_local_num_cells(local_num_cells)
      call set_total_num_cells(local_num_cells)

      allocate (topo%global_indices(local_num_cells))
      do i = 1, topo%local_num_cells
        topo%global_indices(i) = int(topo%graph_conn%vtxdist(irank + 1), ccs_int) + (i - 1)
      end do

      allocate (topo%num_nb(local_num_cells))
      topo%num_nb(:) = max_faces
    end associate

    call nullify_topo_object()

  end subroutine read_topology_topo

  subroutine read_topology_connectivity(shared_env, reader_env, geo_reader, &
                                        global_num_faces, max_faces, &
                                        face_cell1, face_cell1_window, face_cell2, face_cell2_window)

    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: reader_env !< The reader parallel environment
    class(io_process), intent(in) :: geo_reader                                            !< The IO process for reading the file
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

  !v Read the geometry data from an input (HDF5) file
  subroutine read_geometry(shared_env, reader_env, geo_reader, mesh)

    use kinds, only: CCS_MPI_PRECISION

    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: reader_env !< The reader parallel environment
    class(io_process), intent(in) :: geo_reader                                         !< The IO process for reading the file
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
    integer(ccs_int) :: sum_local_num_cells
    integer(ccs_int) :: sum_total_num_cells
    integer(ccs_int) :: max_faces
    integer(ccs_int) :: all_max_faces
    type(face_locator) :: loc_f ! Face locator object
    type(vert_locator) :: loc_v ! Vertex locator object

    integer(ccs_err) :: ierr
    integer :: shared_comm

    integer :: temp_a_f_window, temp_n_f_window, temp_window, temp_x_f_window, temp_x_v_window

    call set_mesh_object(mesh)
    select type (shared_env)
    type is (parallel_environment_mpi)
      shared_comm = shared_env%comm
    class default
      shared_comm = -42
      call error_abort("Unsupported shared environment")
    end select

    call get_max_faces(max_faces)
    if (max_faces == 6) then ! if cell are hexes
      call set_vert_per_cell(8) ! 8 vertices per cell
    else if (max_faces == 4) then ! if cell are tetrahedral
      call set_vert_per_cell(4) ! 4 vertices per cell
    else
      call error_abort("Currently only supporting hex cells.")
    end if

    call get_local_num_cells(local_num_cells)
    call get_total_num_cells(total_num_cells)
    select type (shared_env)
    type is (parallel_environment_mpi)
      call mpi_allreduce(local_num_cells, sum_local_num_cells, 1, MPI_INTEGER, MPI_SUM, shared_env%comm, ierr)
      call mpi_allreduce(total_num_cells, sum_total_num_cells, 1, MPI_INTEGER, MPI_SUM, shared_env%comm, ierr)
      call mpi_allreduce(max_faces, all_max_faces, 1, MPI_INTEGER, MPI_MAX, shared_env%comm, ierr)
    class default
      call error_abort("invalid parallel environment")
    end select

    call get_vert_per_cell(vert_per_cell)

    ! Read attribute "scalefactor"
    if (is_valid(reader_env)) then
      call read_scalar(geo_reader, "scalefactor", mesh%geo%scalefactor)
    end if
    call MPI_Bcast(mesh%geo%scalefactor, 1, CCS_MPI_PRECISION, 0, shared_comm, ierr)

    ! Starting point for reading chunk of data
    vol_p_start = 0

    ! How many data points will be read?
    call get_global_num_cells(global_num_cells)
    vol_p_count = global_num_cells

    associate (local_offset => mesh%topo%shared_array_local_offset, &
               total_offset => mesh%topo%shared_array_total_offset)
      ! Allocate shared memory array for cell volumes
      call create_shared_array(shared_env, sum_total_num_cells, mesh%geo%volumes, mesh%geo%volumes_window)

      ! Read variable "/cell/vol"
      call create_shared_array(shared_env, global_num_cells, temp_vol_c, &
                               temp_window)
      if (is_valid(reader_env)) then
        call read_array(geo_reader, "/cell/vol", vol_p_start, vol_p_count, temp_vol_c)
      end if
      call sync(shared_env)
      mesh%geo%volumes(total_offset + 1:total_offset + total_num_cells) = temp_vol_c(mesh%topo%natural_indices(:))
      call sync(shared_env)
      call destroy_shared_array(shared_env, temp_vol_c, temp_window)

      ! Starting point for reading chunk of data
      x_p_start = [0, 0]

      ! How many data points will be read?
      x_p_count = [ndim, global_num_cells]

      ! Allocate shared memory array for cell centre coordinates
      call create_shared_array(shared_env, [ndim, sum_total_num_cells], mesh%geo%x_p, mesh%geo%x_p_window)

      ! Read variable "/cell/x"
      call create_shared_array(shared_env, [ndim, global_num_cells], temp_x_p, &
                               temp_window)
      if (is_valid(reader_env)) then
        call read_array(geo_reader, "/cell/x", x_p_start, x_p_count, temp_x_p)
      end if
      call sync(shared_env)
      mesh%geo%x_p(:, total_offset + 1:total_offset + total_num_cells) = temp_x_p(:, mesh%topo%natural_indices(:))
      call sync(shared_env)
      call destroy_shared_array(shared_env, temp_x_p, temp_window)

      ! Allocate temporary arrays for face centres, face normals, face areas and vertex coords
      call get_global_num_faces(global_num_faces)

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
      call get_global_num_vertices(global_num_vertices)
      call create_shared_array(shared_env, [ndim, global_num_vertices], temp_x_v, &
                               temp_x_v_window)
      f_xn_count(1) = ndim
      f_xn_count(2) = global_num_vertices
      if (is_valid(reader_env)) then
        call read_array(geo_reader, "/vert", f_xn_start, f_xn_count, temp_x_v)
      end if
      call sync(shared_env)

      ! Allocate shared memory arrays for face centres, face normals, face areas and vertex coordinates
      call create_shared_array(shared_env, [ndim, all_max_faces, sum_local_num_cells], mesh%geo%x_f, mesh%geo%x_f_window)
      call create_shared_array(shared_env, [ndim, all_max_faces, sum_local_num_cells], mesh%geo%face_normals, mesh%geo%face_normals_window) 
      call create_shared_array(shared_env, [all_max_faces, sum_local_num_cells], mesh%geo%face_areas, mesh%geo%face_areas_window)
      call create_shared_array(shared_env, [ndim, vert_per_cell, sum_local_num_cells], mesh%geo%vert_coords, mesh%geo%vert_coords_window)

      ! Procs fill local data
      call sync(shared_env)
      do local_icell = 1, local_num_cells ! loop over cells owned by current process
        call create_cell_locator(local_icell, loc_p)
        call get_natural_index(loc_p, global_icell)

        do j = 1, max_faces ! loop over all faces for each cell
          call create_face_locator(local_icell, j, loc_f)

          n = mesh%topo%global_face_indices(j, global_icell)
          call set_centre(loc_f, temp_x_f(:, n))
          do i = 1, ndim ! loop over dimensions
            ! Map from temp array to mesh for face centres and face normals
            mesh%geo%face_normals(i, j, local_icell + local_offset) = temp_n_f(i, n)
          end do

          ! Map from temp array to mesh for face areas
          call set_area(temp_a_f(n), loc_f)
        end do

        do j = 1, vert_per_cell ! loop over all vertices for each cell
          call create_vert_locator(local_icell, j, loc_v)

          n = mesh%topo%loc_global_vertex_indices(j, local_icell)
          call set_centre(loc_v, temp_x_v(:, n))
        end do

      end do

      ! Correct normal orientations and norms
      do index_p = 1, local_num_cells ! loop over cells owned by current process

        call create_cell_locator(index_p, loc_p)
        call get_centre(loc_p, x_p)
        call count_neighbours(loc_p, nnb)

        do j = 1, nnb ! loop over all faces for each cell

          call create_face_locator(index_p, j, loc_f)
          call get_face_normal(loc_f, face_normal)
          call get_centre(loc_f, x_f)

          if (dot_product(face_normal(:), x_f - x_p) < 0.0_ccs_real) then
            face_normal = -face_normal
          end if

          ! Normalise face normals too
          call set_normal(loc_f, face_normal / norm2(face_normal))

        end do
      end do
    end associate

    ! Delete temp arrays
    call sync(shared_env)
    call destroy_shared_array(shared_env, temp_x_f, temp_x_f_window)
    call destroy_shared_array(shared_env, temp_x_v, temp_x_v_window)
    call destroy_shared_array(shared_env, temp_n_f, temp_n_f_window)
    call destroy_shared_array(shared_env, temp_a_f, temp_a_f_window)

    call compute_face_interpolation(mesh)
    call nullify_mesh_object()

  end subroutine read_geometry

  !v Utility constructor to build a 2D mesh with hex cells.
  !
  !  Builds a Cartesian grid of nx*ny cells.
  function build_square_mesh(par_env, shared_env, run_options) result(mesh)

    use partitioning, only: compute_partitioner_input

    class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared memory environment
    type(ccs_options), intent(in) :: run_options

    character(len=128), dimension(4) :: bnd_names      !< Boundary name list

    type(ccs_mesh) :: mesh                             !< The resulting mesh.

    character(:), allocatable :: error_message

    real(ccs_real) :: side_length
    integer :: cps

    side_length = run_options%mesh%domain_size
    cps = run_options%mesh%cps

    bnd_names = run_options%mesh%bnd_names

    if (cps * cps < par_env%num_procs) then
      error_message = "ERROR: Global number of cells < number of ranks. &
                      &Increase the mesh size or reduce the number of MPI ranks."
      call error_abort(error_message)
    end if

    call set_mesh_object(mesh)
    call set_mesh_generated(.true.)
    call nullify_mesh_object()

    call build_square_topology(par_env, shared_env, cps, mesh)

    call compute_partitioner_input(par_env, shared_env, mesh)

    call mesh_partition_reorder(par_env, shared_env, run_options, mesh)

    call set_offsets(shared_env, mesh)

    call build_square_geometry(par_env, shared_env, cps, side_length, mesh)

    call cleanup_topo(shared_env, mesh)

    ! Create boundary names list
    mesh%bnd_names = bnd_names
    call check_mesh_bnd_names(par_env, mesh)

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
    integer(ccs_int) :: face_index_counter    ! global face counter

    integer(ccs_int) :: global_index_nb ! The global index of a neighbour cell

    integer(ccs_int) :: nglobal          ! The global number of cells
    integer(ccs_int) :: local_num_cells  ! The local number of cells
    integer(ccs_int) :: total_num_cells  ! The total number of cells
    integer(ccs_int) :: global_num_faces ! The global number of faces
    integer(ccs_int) :: max_faces        ! The maximum number of faces per cell
    integer(ccs_int) :: vert_per_cell    ! The number of vertices per cell
    integer(ccs_int) :: global_num_vertices
    integer(ccs_int), dimension(2) :: nb_direction  ! Array indicating direction of neighbour

    type(face_locator) :: loc_f

    integer(ccs_int), dimension(2) :: length

    integer(ccs_int), dimension(:), allocatable :: new_halos ! New halos generated during mesh build

    call set_mesh_object(mesh)
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
        call set_global_num_cells(nglobal)
        call set_global_num_vertices((cps + 1)**2)

        call set_global_num_faces(global_num_faces)
        call set_max_faces(max_faces)

        ! Just to make sure we are working with the same numbers as the mesh object.
        call get_global_num_cells(nglobal)
        call get_global_num_vertices(global_num_vertices)
        call get_global_num_faces(global_num_faces)
        call get_max_faces(max_faces)

        ! Associate aliases to make code easier to read
        associate (h => mesh%geo%h)

          ! Determine ownership range
          start_global = global_start(nglobal, par_env%proc_id, par_env%num_procs)
          local_num_cells = local_count(nglobal, par_env%proc_id, par_env%num_procs)
          call set_local_num_cells(local_num_cells)
          call get_local_num_cells(local_num_cells) ! Ensure using correct value

          ! Abort the execution if any rank has 0 local cells
          if (local_num_cells <= 0) then
            call error_abort("ERROR: Zero local cells found.")
          end if

          call set_total_num_cells(local_num_cells) ! Set initial value
          end_global = start_global + (local_num_cells - 1)

          ! Set number of vertices per cell
          call set_vert_per_cell(4_ccs_int)

          call get_vert_per_cell(vert_per_cell)

          ! Allocate mesh topolgy arrays
          allocate (mesh%topo%global_indices(local_num_cells))
          allocate (mesh%topo%num_nb(local_num_cells))
          allocate (mesh%topo%nb_indices(max_faces, local_num_cells))
          allocate (mesh%topo%face_indices(max_faces, local_num_cells))

          ! Initialise mesh arrays
          mesh%topo%num_nb(:) = max_faces ! All cells have 4 neighbours (possibly ghost/boundary cells)

          ! Initialise neighbour indices
          mesh%topo%nb_indices(:, :) = 0_ccs_int

          ! First set the global index of local cells
          index_counter = 1_ccs_int
          do i = start_global, end_global
            mesh%topo%global_indices(index_counter) = i
            index_counter = index_counter + 1
          end do

          ! Assemble cells and faces
          ! @note Negative neighbour indices are used to indicate boundaries using the same
          !       numbering as cell-relative neighbour indexing, i.e.
          !        -1 = left boundary
          !        -2 = right boundary
          !        -3 = bottom boundary
          !        -4 = top boundary
          index_counter = 1_ccs_int ! Set local indexing starting from 1...n
          allocate (new_halos(0))
          do i = start_global, end_global
            ii = i - 1_ccs_int
            nb_direction = 0_ccs_int

            ! Construct left (1) face/neighbour
            nb_direction(1) = left
            face_counter = left
            call add_neighbour(i, face_counter, index_counter, nb_direction, cps, cps, cps, mesh, new_halos)

            ! Construct right (2) face/neighbour
            nb_direction(1) = right
            face_counter = right
            call add_neighbour(i, face_counter, index_counter, nb_direction, cps, cps, cps, mesh, new_halos)

            ! Construct bottom (3) face/neighbour
            nb_direction(1) = bottom
            face_counter = bottom
            call add_neighbour(i, face_counter, index_counter, nb_direction, cps, cps, cps, mesh, new_halos)

            ! Construct top (4) face/neighbour
            nb_direction(1) = top
            face_counter = top
            call add_neighbour(i, face_counter, index_counter, nb_direction, cps, cps, cps, mesh, new_halos)

            index_counter = index_counter + 1_ccs_int
          end do
        end associate

        ! Append new halos to global indices
        mesh%topo%global_indices = [mesh%topo%global_indices, new_halos]
        deallocate (new_halos)

        call set_total_num_cells(size(mesh%topo%global_indices))
        call get_total_num_cells(total_num_cells)
        call set_halo_num_cells(total_num_cells - local_num_cells)

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

            call create_face_locator(i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
            face_index_counter = face_index_counter + 1_ccs_int
          else
            ! If internal left face, nothing to be done, the face will be linked as a right face from its neighbour
            ! global_index_nb = i - 1_ccs_int
            ! mesh%topo%face_cell1(face_index_counter) = global_index_nb
            ! mesh%topo%face_cell2(face_index_counter) = i
          end if

          ! Construct right (2) face/neighbour
          face_counter = right
          if (modulo(ii, cps) == (cps - 1_ccs_int)) then
            global_index_nb = -right

            call create_face_locator(i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
          else
            global_index_nb = i + 1_ccs_int

            call create_face_locator(i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            call create_face_locator(global_index_nb, left, loc_f)
            call set_global_index(face_index_counter, loc_f)
          end if
          face_index_counter = face_index_counter + 1_ccs_int

          ! Construct bottom (3) face/neighbour
          face_counter = bottom
          if (modulo(ii / cps, cps) == 0_ccs_int) then
            global_index_nb = -bottom

            call create_face_locator(i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
            face_index_counter = face_index_counter + 1_ccs_int
          else
            ! If internal bottom face, nothing to be done, the face will be linked as a top face from its neighbour
            ! global_index_nb = i - nx
            ! mesh%topo%face_cell1(face_index_counter) = global_index_nb
            ! mesh%topo%face_cell2(face_index_counter) = i
          end if

          ! Construct top (4) face/neighbour
          face_counter = top
          if (modulo(ii / cps, cps) == (cps - 1_ccs_int)) then
            global_index_nb = -top

            call create_face_locator(i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
          else
            global_index_nb = i + cps

            call create_face_locator(i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            call create_face_locator(global_index_nb, bottom, loc_f)
            call set_global_index(face_index_counter, loc_f)
          end if
          face_index_counter = face_index_counter + 1_ccs_int

        end do

        call set_num_faces(count_mesh_faces())

        call set_cell_face_indices()

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
    call nullify_mesh_object()

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

  subroutine build_square_geometry(par_env, shared_env, cps, side_length, mesh)

    class(parallel_environment), intent(in) :: par_env !< The parallel environment to construct the mesh.
    class(parallel_environment), intent(in) :: shared_env !< The shared memory environment
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
    integer(ccs_int) :: sum_local_num_cells
    integer(ccs_int) :: sum_total_num_cells
    integer(ccs_int) :: all_max_faces
    integer(ccs_int) :: vert_per_cell

    logical :: is_boundary

    integer(ccs_int) :: index_nb        ! The local index of a neighbour cell

    real(ccs_real), dimension(2) :: x_p ! Cell centre array
    type(cell_locator) :: loc_p         ! Cell locator object

    real(ccs_real), dimension(3) :: x_nb_3 ! Cell centre array of neighbour cell
    real(ccs_real), dimension(2) :: x_nb   ! Cell centre array of neighbour cell
    type(neighbour_locator) :: loc_nb      ! the neighbour locator object.

    real(ccs_real), dimension(2) :: x_f    ! Face centre array
    real(ccs_real), dimension(2) :: normal ! Face normal array
    type(face_locator) :: loc_f            ! Face locator object

    real(ccs_real), dimension(2) :: x_v ! Vertex centre array
    type(vert_locator) :: loc_v         ! Vertex locator object

    integer :: ierr

    call set_mesh_object(mesh)
    select type (par_env)
    type is (parallel_environment_mpi)

      call get_local_num_cells(local_num_cells)
      call get_vert_per_cell(vert_per_cell)

      call get_total_num_cells(total_num_cells)
      call get_max_faces(max_faces)

      select type (shared_env)
      type is (parallel_environment_mpi)
        call mpi_allreduce(local_num_cells, sum_local_num_cells, 1, MPI_INTEGER, MPI_SUM, shared_env%comm, ierr)
        call mpi_allreduce(total_num_cells, sum_total_num_cells, 1, MPI_INTEGER, MPI_SUM, shared_env%comm, ierr)
        call mpi_allreduce(max_faces, all_max_faces, 1, MPI_INTEGER, MPI_MAX, shared_env%comm, ierr)
      class default
        call error_abort("invalid parallel environment")
      end select
      call create_shared_array(shared_env, [ndim, sum_total_num_cells], mesh%geo%x_p, mesh%geo%x_p_window)
      call create_shared_array(shared_env, [ndim, all_max_faces, sum_local_num_cells], mesh%geo%x_f, mesh%geo%x_f_window) !< @note Currently hardcoded as a 2D mesh. @endnote
      call create_shared_array(shared_env, [ndim, all_max_faces, sum_local_num_cells], mesh%geo%face_normals, mesh%geo%face_normals_window) ! Currently hardcoded as a 2D mesh.
      call create_shared_array(shared_env, sum_total_num_cells, mesh%geo%volumes, mesh%geo%volumes_window)
      call create_shared_array(shared_env, [all_max_faces, sum_local_num_cells], mesh%geo%face_areas, mesh%geo%face_areas_window)
      call create_shared_array(shared_env, [ndim, vert_per_cell, sum_local_num_cells], mesh%geo%vert_coords, mesh%geo%vert_coords_window)

      mesh%geo%h = side_length / real(cps, ccs_real)
      if (is_root(shared_env)) then
        mesh%geo%volumes(:) = mesh%geo%h**2 !< @note Mesh is square and 2D @endnote
        mesh%geo%face_normals(:, :, :) = 0.0_ccs_real
        mesh%geo%x_p(:, :) = 0.0_ccs_real
        mesh%geo%x_f(:, :, :) = 0.0_ccs_real
        mesh%geo%face_areas(:, :) = mesh%geo%h  ! Mesh is square and 2D
        mesh%geo%vert_coords(:, :, :) = 0.0_ccs_real
      end if
      call sync(shared_env)

      ! Set cell centre
      associate (h => mesh%geo%h)
        do i = 1_ccs_int, total_num_cells
          call create_cell_locator(i, loc_p)
          call get_natural_index(loc_p, ii)

          x_p(1) = (modulo(ii - 1, cps) + 0.5_ccs_real) * h
          x_p(2) = ((ii - 1) / cps + 0.5_ccs_real) * h

          call set_centre(loc_p, x_p)
        end do
        call sync(shared_env)

        do i = 1_ccs_int, local_num_cells
          call create_cell_locator(i, loc_p)
          call get_centre(loc_p, x_p)

          do face_counter = 1_ccs_int, max_faces

            call create_neighbour_locator(loc_p, face_counter, loc_nb)
            call get_boundary_status(loc_nb, is_boundary)
            call create_face_locator(i, face_counter, loc_f)

            if (.not. is_boundary) then
              ! faces are midway between cell centre and nb cell centre
              call get_centre(loc_nb, x_nb_3)
              x_nb(:) = x_nb_3(1:2) ! hacky fix for issue with resolving get_neighbour_centre in 2D

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
        call sync(shared_env)

        do i = 1_ccs_int, local_num_cells
          call create_cell_locator(i, loc_p)
          call get_centre(loc_p, x_p)

          vertex_counter = front_bottom_left
          call create_vert_locator(i, vertex_counter, loc_v)
          x_v(1) = x_p(1) - 0.5_ccs_real * h
          x_v(2) = x_p(2) - 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = front_bottom_right
          call create_vert_locator(i, vertex_counter, loc_v)
          x_v(1) = x_p(1) + 0.5_ccs_real * h
          x_v(2) = x_p(2) - 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = front_top_left
          call create_vert_locator(i, vertex_counter, loc_v)
          x_v(1) = x_p(1) - 0.5_ccs_real * h
          x_v(2) = x_p(2) + 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = front_top_right
          call create_vert_locator(i, vertex_counter, loc_v)
          x_v(1) = x_p(1) + 0.5_ccs_real * h
          x_v(2) = x_p(2) + 0.5_ccs_real * h
          call set_centre(loc_v, x_v)
        end do
        call sync(shared_env)
      end associate

      call compute_face_interpolation(mesh)

    class default
      call error_abort("Unknown parallel environment type.")

    end select

    call nullify_mesh_object()

  end subroutine build_square_geometry

  !v Utility constructor to build a 3D mesh with hex cells.
  !
  !  Builds a Cartesian grid of nx*ny*nz cells.
  function build_mesh(par_env, shared_env, run_options) result(mesh)

    use partitioning, only: compute_partitioner_input
    use profiler

    class(parallel_environment), allocatable, target, intent(in) :: par_env    !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The shared memory environment
    type(ccs_options), intent(in) :: run_options

    character(len=128), dimension(6) :: bnd_names

    type(ccs_mesh) :: mesh                             !< The resulting mesh.

    character(:), allocatable :: error_message

    real(ccs_real) :: side_length
    integer :: nx, ny, nz

    side_length = run_options%mesh%domain_size
    nx = run_options%mesh%cps
    ny = run_options%mesh%cps
    nz = run_options%mesh%cps

    bnd_names = run_options%mesh%bnd_names

    call set_mesh_object(mesh)
    call set_mesh_generated(.true.)
    call nullify_mesh_object()

    if (.not. (nx .eq. ny .and. ny .eq. nz)) then !< @note Must be a cube (for now) @endnote
      error_message = "Only supporting cubes for now - nx, ny and nz must be the same"
      call error_abort(error_message)
    end if

    if (nx * ny * nz < par_env%num_procs) then
      error_message = "ERROR: Global number of cells < number of ranks. &
                      &Increase the mesh size or reduce the number of MPI ranks."
      call error_abort(error_message)
    end if

    call profiler_begin_region("Build mesh topology")
    call build_topology(par_env, shared_env, nx, ny, nz, mesh)
    call profiler_end_region("Build mesh topology")

    call profiler_begin_region("Compute partitioner input")
    call compute_partitioner_input(par_env, shared_env, mesh)
    call profiler_end_region("Compute partitioner input")

    call mesh_partition_reorder(par_env, shared_env, run_options, mesh)

    call set_offsets(shared_env, mesh)

    call profiler_begin_region("Build mesh geometry")
    call build_geometry(par_env, shared_env, nx, ny, nz, side_length, mesh)
    call profiler_end_region("Build mesh geometry")

    call cleanup_topo(shared_env, mesh)

    ! Create boundary names list
    mesh%bnd_names = bnd_names
    call check_mesh_bnd_names(par_env, mesh)

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
    integer(ccs_int) :: a, b, c, d, e       ! Temporary variables

    integer(ccs_int) :: global_index_nb ! The global index of a neighbour cell
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: total_num_cells
    integer(ccs_int) :: global_num_faces
    integer(ccs_int) :: max_faces
    integer(ccs_int) :: vert_per_cell
    integer(ccs_int) :: global_num_vertices
    integer(ccs_int) :: nglobal
    integer(ccs_int), dimension(3) :: nb_direction  ! Array indicating direction of neighbour

    type(face_locator) :: loc_f

    integer(ccs_int), dimension(2) :: length

    integer(ccs_int), dimension(:), allocatable :: new_halos ! New halos generated during mesh build

    call set_mesh_object(mesh)

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
        call set_global_num_cells(nglobal)
        call set_global_num_vertices((nx + 1) * (ny + 1) * (nz + 1))
        call set_global_num_faces(global_num_faces)
        call set_max_faces(max_faces)

        ! Just to make sure we are working with the same numbers as the mesh object.
        call get_global_num_cells(nglobal)
        call get_global_num_vertices(global_num_vertices)
        call get_global_num_faces(global_num_faces)
        call get_max_faces(max_faces)

        ! Determine ownership range
        start_global = global_start(nglobal, par_env%proc_id, par_env%num_procs)
        local_num_cells = local_count(nglobal, par_env%proc_id, par_env%num_procs)
        call set_local_num_cells(local_num_cells)
        call get_local_num_cells(local_num_cells) ! Ensure using correct value

        ! Abort the execution if any rank has 0 local cells
        if (local_num_cells <= 0) then
          call error_abort("ERROR: Zero local cells found.")
        end if

        call set_total_num_cells(local_num_cells) ! Setting initial value
        end_global = start_global + (local_num_cells - 1)

        ! Set max number of faces (constant, 6)
        call set_max_faces(6_ccs_int)

        ! Set number of vertices per cell (constant, 8)
        call set_vert_per_cell(8)

        call get_max_faces(max_faces)
        call get_vert_per_cell(vert_per_cell)

        ! Allocate mesh arrays
        allocate (mesh%topo%global_indices(local_num_cells))
        allocate (mesh%topo%num_nb(local_num_cells))
        allocate (mesh%topo%nb_indices(max_faces, local_num_cells))
        allocate (mesh%topo%face_indices(max_faces, local_num_cells))

        ! Initialise mesh arrays
        mesh%topo%num_nb(:) = max_faces ! All cells have 6 neighbours (possibly ghost/boundary cells)

        ! Initalise neighbour indices
        mesh%topo%nb_indices(:, :) = 0_ccs_int

        ! First set the global index of local cells
        index_counter = 1_ccs_int
        do i = start_global, end_global
          mesh%topo%global_indices(index_counter) = i
          index_counter = index_counter + 1
        end do

        ! Assemble cells and faces
        ! @note Negative neighbour indices are used to indicate boundaries using the same numbering
        !       as cell-relative neighbour indexing, i.e.
        !        -1 = left boundary
        !        -2 = right boundary
        !        -3 = bottom boundary
        !        -4 = top boundary
        !        -5 = back_boundary
        !        -6 = front_boundary
        index_counter = 1_ccs_int ! Set local indexing starting from 1...n
        allocate (new_halos(0))
        do i = start_global, end_global

          ii = i - 1_ccs_int
          nb_direction(:) = 0_ccs_int

          ! Construct left (1) face/neighbour
          nb_direction(1) = left
          face_counter = left
          call add_neighbour(i, face_counter, index_counter, nb_direction, nx, ny, nz, mesh, new_halos)

          ! Construct right (2) face/neighbour
          nb_direction(1) = right
          face_counter = right
          call add_neighbour(i, face_counter, index_counter, nb_direction, nx, ny, nz, mesh, new_halos)

          ! Construct bottom (3) face/neighbour
          nb_direction(1) = bottom
          face_counter = bottom
          call add_neighbour(i, face_counter, index_counter, nb_direction, nx, ny, nz, mesh, new_halos)

          ! Construct top (4) face/neighbour
          nb_direction(1) = top
          face_counter = top
          call add_neighbour(i, face_counter, index_counter, nb_direction, nx, ny, nz, mesh, new_halos)

          ! Construct back (5) face/neighbour
          nb_direction(1) = back
          face_counter = back
          call add_neighbour(i, face_counter, index_counter, nb_direction, nx, ny, nz, mesh, new_halos)

          ! Construct front (6) face/neighbour
          nb_direction(1) = front
          face_counter = front
          call add_neighbour(i, face_counter, index_counter, nb_direction, nx, ny, nz, mesh, new_halos)

          index_counter = index_counter + 1_ccs_int

        end do

        ! Append new halo indices to global indices
        mesh%topo%global_indices = [mesh%topo%global_indices, new_halos]
        deallocate (new_halos)

        call set_total_num_cells(size(mesh%topo%global_indices))
        call get_total_num_cells(total_num_cells)
        call set_halo_num_cells(total_num_cells - local_num_cells)

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

            call create_face_locator(i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
            face_index_counter = face_index_counter + 1_ccs_int
          else
            ! If internal left face, nothing to be done, the face will be linked as a right face from its neighbour
            ! global_index_nb = i - 1_ccs_int
            ! mesh%topo%face_cell1(face_index_counter) = global_index_nb
            ! mesh%topo%face_cell2(face_index_counter) = i
          end if

          ! Construct right (2) face/neighbour
          face_counter = right
          if (modulo(ii, nx) == (nx - 1_ccs_int)) then
            global_index_nb = -right

            call create_face_locator(i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
          else
            global_index_nb = i + 1_ccs_int

            call create_face_locator(i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            call create_face_locator(global_index_nb, left, loc_f)
            call set_global_index(face_index_counter, loc_f)
          end if
          face_index_counter = face_index_counter + 1_ccs_int

          ! Construct bottom (3) face/neighbour
          face_counter = bottom
          if (modulo(ii / nx, ny) == 0_ccs_int) then
            global_index_nb = -bottom

            call create_face_locator(i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
            face_index_counter = face_index_counter + 1_ccs_int
          else
            ! If internal bottom face, nothing to be done, the face will be linked as a top face from its neighbour
            ! global_index_nb = i - nx
            ! mesh%topo%face_cell1(face_index_counter) = global_index_nb
            ! mesh%topo%face_cell2(face_index_counter) = i
          end if

          ! Construct top (4) face/neighbour
          face_counter = top
          if (modulo(ii / nx, ny) == (ny - 1_ccs_int)) then
            global_index_nb = -top

            call create_face_locator(i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
          else
            global_index_nb = i + nx

            call create_face_locator(i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            call create_face_locator(global_index_nb, bottom, loc_f)
            call set_global_index(face_index_counter, loc_f)
          end if
          face_index_counter = face_index_counter + 1_ccs_int

          ! Construct back (5) face/neighbour
          face_counter = back
          if ((ii / (nx * ny)) == 0_ccs_int) then
            global_index_nb = -back

            call create_face_locator(i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
            face_index_counter = face_index_counter + 1_ccs_int
          else
            ! If internal back face, nothing to be done, the face will be linked as a front face from its neighbour
            ! global_index_nb = i - nx * ny
            ! mesh%topo%face_cell1(face_index_counter) = global_index_nb
            ! mesh%topo%face_cell2(face_index_counter) = i
          end if

          ! Construct front (6) face/neighbour
          face_counter = front
          if ((ii / (nx * ny)) == nz - 1_ccs_int) then
            global_index_nb = -front

            call create_face_locator(i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            mesh%topo%bnd_rid(face_index_counter) = global_index_nb
          else
            global_index_nb = i + nx * ny

            call create_face_locator(i, face_counter, loc_f)
            call set_global_index(face_index_counter, loc_f)
            call create_face_locator(global_index_nb, back, loc_f)
            call set_global_index(face_index_counter, loc_f)
          end if
          face_index_counter = face_index_counter + 1_ccs_int

        end do

        call set_num_faces(count_mesh_faces())

        call set_cell_face_indices()

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
    call nullify_mesh_object()

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
  subroutine build_geometry(par_env, shared_env, nx, ny, nz, side_length, mesh)

    class(parallel_environment), intent(in) :: par_env    !< The parallel environment to construct the mesh.
    class(parallel_environment), intent(in) :: shared_env !< The shared parallel environment.
    integer(ccs_int), intent(in) :: nx                    !< Number of cells in the x direction.
    integer(ccs_int), intent(in) :: ny                    !< Number of cells in the y direction.
    integer(ccs_int), intent(in) :: nz                    !< Number of cells in the z direction.
    real(ccs_real), intent(in) :: side_length             !< The length of the side.

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
    integer(ccs_int) :: sum_local_num_cells
    integer(ccs_int) :: sum_total_num_cells
    integer(ccs_int) :: all_max_faces
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

    integer :: ierr

    associate (foo => nz) ! Silence unused dummy argument
    end associate

    call set_mesh_object(mesh)
    call get_local_num_cells(local_num_cells) ! Ensure using correct value
    call get_vert_per_cell(vert_per_cell)

    select type (par_env)
    type is (parallel_environment_mpi)

      call get_total_num_cells(total_num_cells)
      call get_max_faces(max_faces)

      select type (shared_env)
      type is (parallel_environment_mpi)
        call mpi_allreduce(local_num_cells, sum_local_num_cells, 1, MPI_INTEGER, MPI_SUM, shared_env%comm, ierr)
        call mpi_allreduce(total_num_cells, sum_total_num_cells, 1, MPI_INTEGER, MPI_SUM, shared_env%comm, ierr)
        call mpi_allreduce(max_faces, all_max_faces, 1, MPI_INTEGER, MPI_MAX, shared_env%comm, ierr)
      class default
        call error_abort("invalid parallel environment")
      end select
      call create_shared_array(shared_env, [ndim, sum_total_num_cells], mesh%geo%x_p, mesh%geo%x_p_window)
      call create_shared_array(shared_env, [ndim, all_max_faces, sum_local_num_cells], mesh%geo%x_f, mesh%geo%x_f_window) !< @note Currently hardcoded as a 2D mesh. @endnote
      call create_shared_array(shared_env, [ndim, all_max_faces, sum_local_num_cells], mesh%geo%face_normals, mesh%geo%face_normals_window) ! Currently hardcoded as a 2D mesh.
      call create_shared_array(shared_env, sum_total_num_cells, mesh%geo%volumes, mesh%geo%volumes_window)
      call create_shared_array(shared_env, [all_max_faces, sum_local_num_cells], mesh%geo%face_areas, mesh%geo%face_areas_window)
      call create_shared_array(shared_env, [ndim, vert_per_cell, sum_local_num_cells], mesh%geo%vert_coords, mesh%geo%vert_coords_window)

      mesh%geo%h = side_length / real(nx, ccs_real) !< @note Assumes cube @endnote
      if (is_root(shared_env)) then
        mesh%geo%volumes(:) = mesh%geo%h**3 !< @note Mesh is cube @endnote
        mesh%geo%face_normals(:, :, :) = 0.0_ccs_real
        mesh%geo%x_p(:, :) = 0.0_ccs_real
        mesh%geo%x_f(:, :, :) = 0.0_ccs_real
        mesh%geo%face_areas(:, :) = mesh%geo%h**2
        mesh%geo%vert_coords(:, :, :) = 0.0_ccs_real
      end if
      call sync(shared_env)

      ! Set cell centre
      associate (h => mesh%geo%h)
        do i = 1_ccs_int, total_num_cells
          call create_cell_locator(i, loc_p)
          call get_natural_index(loc_p, ii)

          x_p(1) = (modulo(ii - 1, nx) + 0.5_ccs_real) * h
          x_p(2) = (modulo((ii - 1) / nx, ny) + 0.5_ccs_real) * h
          x_p(3) = (((ii - 1) / (nx * ny)) + 0.5_ccs_real) * h

          call set_centre(loc_p, x_p)
        end do
        call sync(shared_env)

        do i = 1_ccs_int, local_num_cells
          call create_cell_locator(i, loc_p)
          call get_centre(loc_p, x_p)

          do face_counter = 1_ccs_int, max_faces

            call create_neighbour_locator(loc_p, face_counter, loc_nb)
            call get_boundary_status(loc_nb, is_boundary)
            call create_face_locator(i, face_counter, loc_f)

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
        call sync(shared_env)

        do i = 1_ccs_int, local_num_cells
          call create_cell_locator(i, loc_p)
          call get_centre(loc_p, x_p)

          vertex_counter = front_bottom_left
          call create_vert_locator(i, vertex_counter, loc_v)
          x_v(1) = x_p(1) - 0.5_ccs_real * h
          x_v(2) = x_p(2) - 0.5_ccs_real * h
          x_v(3) = x_p(3) + 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = front_bottom_right
          call create_vert_locator(i, vertex_counter, loc_v)
          x_v(1) = x_p(1) + 0.5_ccs_real * h
          x_v(2) = x_p(2) - 0.5_ccs_real * h
          x_v(3) = x_p(3) + 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = front_top_left
          call create_vert_locator(i, vertex_counter, loc_v)
          x_v(1) = x_p(1) - 0.5_ccs_real * h
          x_v(2) = x_p(2) + 0.5_ccs_real * h
          x_v(3) = x_p(3) + 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = front_top_right
          call create_vert_locator(i, vertex_counter, loc_v)
          x_v(1) = x_p(1) + 0.5_ccs_real * h
          x_v(2) = x_p(2) + 0.5_ccs_real * h
          x_v(3) = x_p(3) + 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = back_bottom_left
          call create_vert_locator(i, vertex_counter, loc_v)
          x_v(1) = x_p(1) - 0.5_ccs_real * h
          x_v(2) = x_p(2) - 0.5_ccs_real * h
          x_v(3) = x_p(3) - 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = back_bottom_right
          call create_vert_locator(i, vertex_counter, loc_v)
          x_v(1) = x_p(1) + 0.5_ccs_real * h
          x_v(2) = x_p(2) - 0.5_ccs_real * h
          x_v(3) = x_p(3) - 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = back_top_left
          call create_vert_locator(i, vertex_counter, loc_v)
          x_v(1) = x_p(1) - 0.5_ccs_real * h
          x_v(2) = x_p(2) + 0.5_ccs_real * h
          x_v(3) = x_p(3) - 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = back_top_right
          call create_vert_locator(i, vertex_counter, loc_v)
          x_v(1) = x_p(1) + 0.5_ccs_real * h
          x_v(2) = x_p(2) + 0.5_ccs_real * h
          x_v(3) = x_p(3) - 0.5_ccs_real * h
          call set_centre(loc_v, x_v)
        end do
        call sync(shared_env)
      end associate

      call compute_face_interpolation(mesh)

    class default
      call error_abort("Unknown parallel environment type.")

    end select
    call nullify_mesh_object()

  end subroutine build_geometry

  !> Helper subroutine to appropriately set local and global neighbour indices
  subroutine add_neighbour(index_p, nb_counter, index_counter, direction, nx, ny, nz, mesh, new_halos)
    integer(ccs_int), intent(in) :: index_p                   !< Global index of cell whose neighbours we're assembling
    integer(ccs_int), intent(in) :: nb_counter                !< the cell-relative index neighbour index
    integer(ccs_int), intent(in) :: index_counter             !< local index of cell whose neighbours we're assembling
    integer(ccs_int), dimension(:), intent(in) :: direction   !< Array containing the direction of the neighbour
    !< relative to the cell
    integer(ccs_int), intent(in) :: nx                        !< Mesh size in x direction
    integer(ccs_int), intent(in) :: ny                        !< Mesh size in y direction
    integer(ccs_int), intent(in) :: nz                        !< Mesh size in z direction
    type(ccs_mesh), intent(inout) :: mesh                     !< The mesh
    integer(ccs_int), dimension(:), allocatable, intent(inout) :: new_halos !< New halo indices

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

    call build_local_mesh_add_neighbour(index_counter, nb_counter, index_nb, global_index_nb, mesh, new_halos)
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
  subroutine build_local_mesh_add_neighbour(index_p, index_p_nb, index_nb, global_index_nb, mesh, new_halos)

    integer(ccs_int), intent(in) :: index_p !< the index of the cell whose neighbours we are assembling
    integer(ccs_int), intent(in) :: index_p_nb !< the cell-relative neighbour index
    integer(ccs_int), intent(in) :: index_nb !< the local index of the neighbour cell
    integer(ccs_int), intent(in) :: global_index_nb !< the global index of the neighbour cell
    type(ccs_mesh), intent(inout) :: mesh !< the mesh we are assembling neighbours on
    integer(ccs_int), dimension(:), allocatable, intent(inout) :: new_halos !< New halo indices

    integer(ccs_int) :: local_num_cells ! The number of local cells
    integer(ccs_int) :: global_num_cells ! The number of global cells

    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb

    call get_local_num_cells(local_num_cells)
    call get_global_num_cells(global_num_cells)

    call create_cell_locator(index_p, loc_p)
    call create_neighbour_locator(loc_p, index_p_nb, loc_nb)
    if ((index_nb >= 1_ccs_int) .and. (index_nb <= local_num_cells)) then
      call set_local_index(index_nb, loc_nb)
    else if (global_index_nb < 0_ccs_int) then
      ! Boundary "neighbour" - local index should also be -ve
      if (.not. (index_nb < 0_ccs_int)) then
        call error_abort("ERROR: boundary neighbours should have -ve indices.")
      end if
      call set_local_index(index_nb, loc_nb)
    else
      ! Neighbour is in a halo
      call add_halo_neighbour(global_index_nb, loc_nb, mesh, new_halos)
    end if

  end subroutine build_local_mesh_add_neighbour

  !v Helper subroutine to add a neighbour in the set of halo cells.
  !
  !  Given the global index of the neighbour:
  !  1) try to locate it at the end of the global indices array - this is for historical reasons, in
  !     the new design only the global indices of local cells should be stored here at this point.
  !  2) check if we've already found this halo neighbour in the new halos list.
  !  3) if the above fails then it is a new halo neighbour.
  !
  !  The local cells and halo cells are mainained separately at this point as it means we can
  !  (relatively) cheaply append halo cells onto the (shorter) new halos list and only at the end
  !  does an expensive concatenation of local and halo cells occur.
  !
  !  This subroutine should only be called after determining a neighbour is a halo (i.e. non-local)
  !  cell.
  subroutine add_halo_neighbour(global_index_nb, loc_nb, mesh, new_halos)

    integer(ccs_int), intent(in) :: global_index_nb !< the global index of the neighbour cell
    type(neighbour_locator), intent(inout) :: loc_nb
    type(ccs_mesh), intent(inout) :: mesh !< the mesh we are assembling neighbours on
    integer(ccs_int), dimension(:), allocatable, intent(inout) :: new_halos !< New halo indices

    integer(ccs_int) :: ng  ! The current number of cells (total = local + halos)
    logical :: found        ! Indicates whether a halo cell was already present
    integer(ccs_int) :: i   ! Cell iteration counter

    integer(ccs_int) :: local_num_cells ! The number of local cells
    integer(ccs_int) :: global_num_cells ! The number of global cells

    call get_local_num_cells(local_num_cells)
    call get_global_num_cells(global_num_cells)

    ! First check if neighbour is already present in halo
    ng = size(mesh%topo%global_indices)
    found = .false.

    ! First look within existing mesh halos
    if (ng > local_num_cells) then
      i = findloc(mesh%topo%global_indices(local_num_cells + 1:ng), global_index_nb, dim=1)
    else
      i = 0
    end if

    if (i > 0) then ! Found neighbour in halos
      found = .true.

      i = i + local_num_cells ! Offset
      call set_local_index(i, loc_nb)
    end if

    if (.not. found) then
      ! Have we seen this halo before?

      ! XXX: abstract
      i = findloc(new_halos, global_index_nb, dim=1)
      if (i > 0) then ! Found neighbour in halos
        found = .true.

        i = i + ng ! Offset
        call set_local_index(i, loc_nb)
      end if
    end if

    if (.not. found) then
      ! Halo is unseen
      call add_new_halo_neighbour(global_index_nb, loc_nb, mesh, new_halos)
    end if

  end subroutine add_halo_neighbour

  ! If neighbour was not present append to global index list (the end of the global index list
  ! becoming its local index).
  ! XXX: Note this currently copies into an n+1 temporary, reallocates and then copies back to
  !      the (extended) original array.
  subroutine add_new_halo_neighbour(global_index_nb, loc_nb, mesh, new_halos)

    integer(ccs_int), intent(in) :: global_index_nb !< the global index of the neighbour cell
    type(neighbour_locator), intent(inout) :: loc_nb
    type(ccs_mesh), intent(inout) :: mesh !< the mesh we are assembling neighbours on
    integer(ccs_int), dimension(:), allocatable, intent(inout) :: new_halos !< New halo indices

    integer(ccs_int) :: ng  ! The current number of cells (total = local + halos)

    integer(ccs_int) :: global_num_cells ! The number of global cells
    integer(ccs_int) :: total_num_cells

    call get_global_num_cells(global_num_cells)

    ng = size(mesh%topo%global_indices) + size(new_halos)

    if ((ng + 1) > global_num_cells) then
      call error_abort("ERROR: Trying to create halo that exceeds global mesh size.")
    end if

    ng = ng + 1
    call set_local_index(ng, loc_nb)
    new_halos = [new_halos, global_index_nb]

    ! Increment total cell count
    call get_total_num_cells(total_num_cells)
    call set_total_num_cells(total_num_cells + 1)

    call get_total_num_cells(total_num_cells)
    if (total_num_cells /= (size(mesh%topo%global_indices) + size(new_halos))) then
      write (log_unit_out, *) total_num_cells, size(mesh%topo%global_indices), size(new_halos)
      call error_abort("ERROR: Local total cell count and size of global indices + new halos not in agreement")
    end if

  end subroutine add_new_halo_neighbour

  !v Count the number of faces in the mesh
  function count_mesh_faces() result(nfaces)

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
    call get_local_num_cells(local_num_cells)
    do index_p = 1, local_num_cells
      call create_cell_locator(index_p, loc_p)
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
  subroutine set_cell_face_indices()

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
    call get_local_num_cells(local_num_cells)
    do index_p = 1, local_num_cells
      call create_cell_locator(index_p, loc_p)
      call count_neighbours(loc_p, nnb)

      do j = 1, nnb
        call create_neighbour_locator(loc_p, j, loc_nb)
        call get_local_index(loc_nb, index_nb)
        call get_boundary_status(loc_nb, is_boundary)

        if (.not. is_boundary) then
          ! Cell with lowest local index assigns an index to the face
          if (index_p < index_nb) then
            face_counter = face_counter + 1
            call set_face_index(index_p, j, face_counter)
          else
            ! Find corresponding face in neighbour cell
            ! (To be improved, this seems inefficient)
            index_f = get_neighbour_face_index(index_p, index_nb)
            call set_face_index(index_p, j, index_f)
          end if
        else
          face_counter = face_counter + 1
          call set_face_index(index_p, j, face_counter)
        end if
      end do  ! End loop over current cell's neighbours
    end do    ! End loop over local cells

  end subroutine set_cell_face_indices

  !v Computes the index of the face shared by the cells denoted by the specified
  !  local index and neighbouring index
  function get_neighbour_face_index(index_p, index_nb) result(index_f)
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

    call create_cell_locator(index_nb, loc_nb)
    call count_neighbours(loc_nb, nnb_nb)
    do k = 1, nnb_nb
      call create_neighbour_locator(loc_nb, k, loc_nb_nb)
      call get_local_index(loc_nb_nb, index_nb_nb)
      if (index_nb_nb == index_p) then
        call create_face_locator(index_nb, k, loc_f)
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
    real(ccs_real), dimension(ndim) :: normal ! face normal

    integer(ccs_int) :: index_f
    integer(ccs_int) :: global_p, global_nb
    integer(ccs_int) :: natural_p, natural_nb
    integer(ccs_int) :: global_f
    integer :: rank, ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    
    if (allocated(mesh%geo%face_interpol)) then
      deallocate (mesh%geo%face_interpol)
    end if

    call get_num_faces(num_faces)
    allocate (mesh%geo%face_interpol(num_faces))

    ! Safe guard to make sure we go through all faces
    mesh%geo%face_interpol(:) = -1.0_ccs_real

    call get_local_num_cells(local_num_cells)

    do index_p = 1, local_num_cells
      call create_cell_locator(index_p, loc_p)
      call count_neighbours(loc_p, nnb)

      do j = 1, nnb
        call create_neighbour_locator(loc_p, j, loc_nb)
        call get_boundary_status(loc_nb, is_boundary)
        call create_face_locator(index_p, j, loc_f)

        if (.not. is_boundary) then
          call get_local_index(loc_nb, index_nb)

          call get_centre(loc_f, x_f)
          call get_centre(loc_nb, x_nb)
          call get_centre(loc_p, x_p)
          call get_face_normal(loc_f, normal)

          ! Project cell centres (P and N) onto face normal (-> P' and N') and compute the ratio |fP'|/ (|fP'|+ |fN'|)
          ! This is equivalent (via Thales' theorem) to getting the ratio between f'P over PN where f' is the point on the face intersecting NP
          interpol_factor = dot_product(normal, x_f - x_p) / abs(dot_product(normal, x_f - x_p) - dot_product(normal, x_f - x_nb))

          ! Face interpolation diagnostics
          if (interpol_factor < 0.0_ccs_real .or. interpol_factor > 1.0_ccs_real) then
            call get_local_index(loc_f, index_f)
            call get_global_index(loc_p, global_p)
            call get_global_index(loc_nb, global_nb)
            call get_natural_index(loc_p, natural_p)
            call get_natural_index(loc_nb, natural_nb)

            global_f = -1_ccs_int
            if (associated(mesh%topo%global_face_indices)) then
              global_f = mesh%topo%global_face_indices(j, natural_p)
            end if

            write (log_unit_out, *) "FACE INTERPOLATION DIAGNOSTIC"
            write (log_unit_out, *) "rank:", rank
            write (log_unit_out, *) "face slot/local/global:", j, index_f, global_f
            write (log_unit_out, *) "cell local/global/natural:", index_p, global_p, natural_p
            write (log_unit_out, *) "neighbour local/global/natural:", index_nb, global_nb, natural_nb
            write (log_unit_out, *) "x_cell:", x_p
            write (log_unit_out, *) "x_neighbour:", x_nb
            write (log_unit_out, *) "x_face:", x_f
            write (log_unit_out, *) "interpolation factor:", interpol_factor
          end if
          ! End face interpolation diagnostics

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

       ! Face interpolation diagnostics
       do index_f = 1, num_faces
          if (mesh%geo%face_interpol(index_f) < 0.0_ccs_real .or. &
               mesh%geo%face_interpol(index_f) > 1.0_ccs_real) then
             write (log_unit_out, *) &
                  "Invalid stored face: rank/local face/value:", &
                  rank, index_f, mesh%geo%face_interpol(index_f)
          end if
       end do
       flush (log_unit_out)
       ! End face interpolation diagnostics

       call error_abort("Face interpolation out of bound.")
    end if

  end subroutine

  ! Populate mesh%topo%graph_conn%global_partition with a split of cells in stride using global_start and local_count
  subroutine partition_stride(par_env, mesh)
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(ccs_mesh), intent(inout) :: mesh                             !< The resulting mesh.

    mesh%topo%graph_conn%local_partition(:) = par_env%proc_id

  end subroutine partition_stride

  !v Compute the global start index of a process if each receives an equal share of the items,
  !  ordered by process ID. Any remainder items are distributed evenly between the remainder lower
  !  process IDs.
  pure integer function global_start(n, procid, nproc)

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
  pure integer function local_count(n, procid, nproc)

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

    write (log_unit_out, *) par_env%proc_id, "############################# Print Geometry ########################################"

    write (log_unit_out, *) par_env%proc_id, "h                  : ", mesh%geo%h
    write (log_unit_out, *) par_env%proc_id, "scalefactor        : ", mesh%geo%scalefactor
    write (log_unit_out, *) ""

    associate (local_offset => mesh%topo%shared_array_local_offset, &
               total_offset => mesh%topo%shared_array_total_offset)

      if (associated(mesh%geo%volumes)) then
        write (log_unit_out, *) par_env%proc_id, "volumes     : ", mesh%geo%volumes(1 + total_offset:nb_elem + total_offset)
      else
        write (log_unit_out, *) par_env%proc_id, "volumes     : UNALLOCATED"
      end if

      if (allocated(mesh%geo%face_interpol)) then
        write (log_unit_out, *) par_env%proc_id, "face_interpol          : ", mesh%geo%face_interpol(1:nb_elem)
      else
        write (log_unit_out, *) par_env%proc_id, "face_interpol          : UNALLOCATED"
      end if

      write (log_unit_out, *) ""
      if (associated(mesh%geo%face_areas)) then
        do i = 1, nb_elem
          write (log_unit_out, *) par_env%proc_id, "face_areas(1:" // str(nb_elem / 2) // ", " // str(i) // ")", &
            mesh%geo%face_areas(1:nb_elem / 2, i + local_offset)
        end do
      else
        write (log_unit_out, *) par_env%proc_id, "face_areas             : UNALLOCATED"
      end if

      write (log_unit_out, *) ""
      if (associated(mesh%geo%x_p)) then
        do i = 1, nb_elem
          write (log_unit_out, *) par_env%proc_id, "x_p(:)", mesh%geo%x_p(:, i + total_offset)
        end do
      else
        write (log_unit_out, *) par_env%proc_id, "x_p                    : UNALLOCATED"
      end if

      write (log_unit_out, *) ""
      if (associated(mesh%geo%x_f)) then
        do i = 1, nb_elem
          write (log_unit_out, *) par_env%proc_id, "x_f(2, 1:" // str(nb_elem / 2) // ", " // str(i) // ")", &
            mesh%geo%x_f(2, 1:nb_elem / 2, i + local_offset)
        end do
      else
        write (log_unit_out, *) par_env%proc_id, "x_f                    : UNALLOCATED"
      end if

      write (log_unit_out, *) ""
      if (associated(mesh%geo%face_normals)) then
        do i = 1, nb_elem
          write (log_unit_out, *) par_env%proc_id, "face_normals(2, 1:" // str(nb_elem / 2) // ", " // str(i) // ")", &
            mesh%geo%face_normals(2, 1:nb_elem / 2, i + local_offset)
        end do
      else
        write (log_unit_out, *) par_env%proc_id, "face_normals          : UNALLOCATED"
      end if

      write (log_unit_out, *) ""
      if (associated(mesh%geo%vert_coords)) then
        do i = 1, nb_elem
          write (log_unit_out, *) par_env%proc_id, "vert_coords(2, 1:" // str(nb_elem / 2) // ", " // str(i) // ")", &
            mesh%geo%vert_coords(2, 1:nb_elem / 2, i + local_offset)
        end do
      else
        write (log_unit_out, *) par_env%proc_id, "vert_coords           : UNALLOCATED"
      end if
    end associate

    write (log_unit_out, *) par_env%proc_id, "############################# End Print Geometry ########################################"

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

    write (log_unit_out, *) par_env%proc_id, "############################# Print Topology ########################################"

    call get_local_num_cells(local_num_cells)
    call get_global_num_cells(global_num_cells)
    call get_halo_num_cells(halo_num_cells)
    call get_total_num_cells(total_num_cells)
    call get_global_num_faces(global_num_faces)
    call get_num_faces(num_faces)
    call get_max_faces(max_faces)
    call get_vert_per_cell(vert_per_cell)
    call get_global_num_vertices(global_num_vertices)

    write (log_unit_out, *) par_env%proc_id, "global_num_cells    : ", global_num_cells
    write (log_unit_out, *) par_env%proc_id, "local_num_cells     : ", local_num_cells
    write (log_unit_out, *) par_env%proc_id, "halo_num_cells      : ", halo_num_cells
    write (log_unit_out, *) par_env%proc_id, "total_num_cells     : ", total_num_cells
    write (log_unit_out, *) par_env%proc_id, "global_num_vertices : ", global_num_vertices
    write (log_unit_out, *) par_env%proc_id, "vert_per_cell       : ", vert_per_cell
    write (log_unit_out, *) par_env%proc_id, "global_num_faces    : ", global_num_faces
    write (log_unit_out, *) par_env%proc_id, "num_faces           : ", num_faces
    write (log_unit_out, *) par_env%proc_id, "max_faces           : ", max_faces
    write (log_unit_out, *) ""

    if (allocated(mesh%topo%global_indices)) then
      write (log_unit_out, *) par_env%proc_id, "global_indices     : ", mesh%topo%global_indices(1:nb_elem)
    else
      write (log_unit_out, *) par_env%proc_id, "global_indices     : UNALLOCATED"
    end if

    if (allocated(mesh%topo%num_nb)) then
      write (log_unit_out, *) par_env%proc_id, "num_nb             : ", mesh%topo%num_nb(1:nb_elem)
    else
      write (log_unit_out, *) par_env%proc_id, "num_nb             : UNALLOCATED"
    end if

    if (associated(mesh%topo%face_cell1)) then
      write (log_unit_out, *) par_env%proc_id, "face_cell1        : ", mesh%topo%face_cell1(1:nb_elem)
    else
      write (log_unit_out, *) par_env%proc_id, "face_cell1        : UNALLOCATED"
    end if

    if (associated(mesh%topo%face_cell2)) then
      write (log_unit_out, *) par_env%proc_id, "face_cell2        : ", mesh%topo%face_cell2(1:nb_elem)
    else
      write (log_unit_out, *) par_env%proc_id, "face_cell2        : UNALLOCATED"
    end if

    if (associated(mesh%topo%bnd_rid)) then
      write (log_unit_out, *) par_env%proc_id, "bnd_rid           : ", mesh%topo%bnd_rid(1:nb_elem)
    else
      write (log_unit_out, *) par_env%proc_id, "bnd_rid           : UNALLOCATED"
    end if

    if (allocated(mesh%topo%graph_conn%vwgt)) then
      write (log_unit_out, *) par_env%proc_id, "vwgt              : ", mesh%topo%graph_conn%vwgt(1:nb_elem)
    else
      write (log_unit_out, *) par_env%proc_id, "vwgt              : UNALLOCATED"
    end if

    if (allocated(mesh%topo%graph_conn%adjwgt)) then
      write (log_unit_out, *) par_env%proc_id, "adjwgt            : ", mesh%topo%graph_conn%adjwgt(1:nb_elem)
    else
      write (log_unit_out, *) par_env%proc_id, "adjwgt            : UNALLOCATED"
    end if

    if (allocated(mesh%topo%graph_conn%local_partition)) then
      write (log_unit_out, *) par_env%proc_id, "local_partition   : ", mesh%topo%graph_conn%local_partition(1:nb_elem)
    else
      write (log_unit_out, *) par_env%proc_id, "local_partition   : UNALLOCATED"
    end if

    write (log_unit_out, *) ""
    if (associated(mesh%topo%global_face_indices)) then
      do i = 1, nb_elem
        write(log_unit_out,*) par_env%proc_id, "global_face_indices(1:"    //    str(nb_elem / 2)    //    ", "    //    str(i)    //    ")", mesh%topo%global_face_indices(1:nb_elem / 2, i)
      end do
    else
      write (log_unit_out, *) par_env%proc_id, "global_face_indices   : UNALLOCATED"
    end if

    write (log_unit_out, *) ""
    if (associated(mesh%topo%global_vertex_indices)) then
      do i = 1, nb_elem
        write(log_unit_out,*) par_env%proc_id, "global_vertex_indices(1:"    //    str(nb_elem / 2)    //    ", "    //    str(i)    //    ")", mesh%topo%global_vertex_indices(1:nb_elem / 2, i)
      end do
    else
      write (log_unit_out, *) par_env%proc_id, "global_vertex_indices : UNALLOCATED"
    end if

    write (log_unit_out, *) ""
    if (allocated(mesh%topo%face_indices)) then
      do i = 1, nb_elem
    write(log_unit_out,*) par_env%proc_id, "face_indices(1:"  //  str(nb_elem / 2)  //  ", "  //  str(i)  //  ")", mesh%topo%face_indices(1:nb_elem / 2, i)
      end do
    else
      write (log_unit_out, *) par_env%proc_id, "face_indices          : UNALLOCATED"
    end if

    write (log_unit_out, *) ""
    if (allocated(mesh%topo%nb_indices)) then
      do i = 1, nb_elem
        write(log_unit_out,*) par_env%proc_id, "nb_indices(1:"  //  str(nb_elem / 2)  //  ", "  //  str(i)  //  ")", mesh%topo%nb_indices(1:nb_elem / 2, i)
      end do
    else
      write (log_unit_out, *) par_env%proc_id, "nb_indices            : UNALLOCATED"
    end if

    write (log_unit_out, *) par_env%proc_id, "############################# End Print Topology ########################################"

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
  end subroutine cleanup_topo

  subroutine destroy_mesh_storage(par_env, mesh, owns_shared_arrays)

    use partitioning, only: cleanup_partitioner_data

    type(ccs_mesh), target, intent(inout) :: mesh   !< The mesh
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    logical, intent(in), optional :: owns_shared_arrays !< Whether this mesh owns its shared arrays

    logical :: owns_shared

    owns_shared = .true.
    if (present(owns_shared_arrays)) owns_shared = owns_shared_arrays

    call cleanup_partitioner_data(mesh)

    if (owns_shared) then
      call cleanup_topo(par_env, mesh)

      if (associated(mesh%geo%face_areas)) then
        call destroy_shared_array(par_env, mesh%geo%face_areas, mesh%geo%face_areas_window)
        call dprint("mesh%geo%face_areas deallocated.")
      end if

      if (associated(mesh%geo%volumes)) then
        call destroy_shared_array(par_env, mesh%geo%volumes, mesh%geo%volumes_window)
        call dprint("mesh%geo%volumes deallocated.")
      end if

      if (associated(mesh%geo%x_p)) then
        call destroy_shared_array(par_env, mesh%geo%x_p, mesh%geo%x_p_window)
        call dprint("mesh%geo%x_p deallocated.")
      end if

      if (associated(mesh%geo%x_f)) then
        call destroy_shared_array(par_env, mesh%geo%x_f, mesh%geo%x_f_window)
        call dprint("mesh%geo%x_f deallocated.")
      end if

      if (associated(mesh%geo%face_normals)) then
        call destroy_shared_array(par_env, mesh%geo%face_normals, mesh%geo%face_normals_window)
        call dprint("mesh%geo%face_normals deallocated.")
      end if

      if (associated(mesh%geo%vert_coords)) then
        call destroy_shared_array(par_env, mesh%geo%vert_coords, mesh%geo%vert_coords_window)
        call dprint("mesh%geo%vert_coords deallocated.")
      end if
    else
      nullify (mesh%topo%global_face_indices)
      nullify (mesh%topo%global_vertex_indices)
      nullify (mesh%topo%face_cell1)
      nullify (mesh%topo%face_cell2)
      nullify (mesh%topo%bnd_rid)
      nullify (mesh%geo%face_areas)
      nullify (mesh%geo%volumes)
      nullify (mesh%geo%x_p)
      nullify (mesh%geo%x_f)
      nullify (mesh%geo%face_normals)
      nullify (mesh%geo%vert_coords)
    end if

    if (allocated(mesh%topo%natural_indices)) then
      deallocate (mesh%topo%natural_indices)
      call dprint("mesh%topo%natural_indices deallocated.")
    end if

    if (allocated(mesh%topo%global_indices)) then
      deallocate (mesh%topo%global_indices)
      call dprint("mesh%topo%global_indices deallocated.")
    end if

    if (allocated(mesh%topo%loc_global_vertex_indices)) then
      deallocate (mesh%topo%loc_global_vertex_indices)
      call dprint("mesh%topo%loc_global_vertex_indices deallocated.")
    end if

    if (allocated(mesh%topo%face_indices)) then
      deallocate (mesh%topo%face_indices)
      call dprint("mesh%topo%face_indices deallocated.")
    end if

    if (allocated(mesh%topo%nb_indices)) then
      deallocate (mesh%topo%nb_indices)
      call dprint("mesh%topo%nb_indices deallocated.")
    end if

    if (allocated(mesh%topo%num_nb)) then
      deallocate (mesh%topo%num_nb)
      call dprint("mesh%topo%num_nb deallocated.")
    end if

    if (allocated(mesh%geo%face_interpol)) then
      deallocate (mesh%geo%face_interpol)
      call dprint("mesh%geo%face_interpol deallocated.")
    end if

    if (allocated(mesh%bnd_names)) then
      deallocate (mesh%bnd_names)
      call dprint("mesh%bnd_names deallocated.")
    end if

    mesh = ccs_mesh()

  end subroutine destroy_mesh_storage

  subroutine mesh_partition_reorder(par_env, shared_env, run_options, mesh)

    use partitioning, only: partition_kway, &
                            compute_connectivity, &
                            compute_partitioner_input, &
                            cleanup_partitioner_data, &
                            print_partition_quality
    use profiler

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    class(parallel_environment), allocatable, target, intent(in) :: shared_env !< The parallel environment
    type(ccs_options), intent(in) :: run_options
    type(ccs_mesh), intent(inout) :: mesh                                   !< The mesh

    call set_mesh_object(mesh)

    if (par_env%num_procs > 1) then
      call profiler_begin_region("Kway partitioning")
      call partition_kway(par_env, mesh)
      call profiler_end_region("Kway partitioning")
    else
      call profiler_begin_region("Strided partitioning")
      call partition_stride(par_env, mesh)
      call profiler_end_region("Strided partitioning")
    end if

    call sync(shared_env)

    call print_partition_quality(par_env, run_options)

    call profiler_begin_region("Computing connectivity")
    call compute_connectivity(par_env, shared_env, mesh)
    call profiler_end_region("Computing connectivity")

! insert halo / local cells computation here

    call print_bandwidth(par_env, run_options)

    call profiler_begin_region("Reordering")
    call reorder_cells(par_env, shared_env, mesh)
    call profiler_end_region("Reordering")

    call cleanup_partitioner_data(mesh)

    call print_bandwidth(par_env, run_options)

    call nullify_mesh_object()

  end subroutine mesh_partition_reorder

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

  ! Build adjacency matrix for local cells
  pure subroutine build_adjacency_matrix(xadj, adjncy)

    integer(ccs_int), allocatable, dimension(:), intent(out) :: xadj   !< Array that points to where in adjncy
    !  the list for each cell begins and ends
    integer(ccs_int), allocatable, dimension(:), intent(out) :: adjncy !< Array storing adjacency lists for each cell consecutively

    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: ctr
    integer(ccs_int) :: idx
    integer(ccs_int) :: i, j, nnb
    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb
    logical :: cell_local

    integer, parameter :: alloc_step = 100 ! How many entries to increase an array by
    integer, dimension(:), allocatable :: tmp

    call get_local_num_cells(local_num_cells)
    allocate (xadj(local_num_cells + 1))
    allocate (adjncy(3 * local_num_cells)) ! A reasonable starting point (minimal faces per cell is 3 == TRI)
    ctr = 1
    xadj(1) = ctr

    do i = 1, local_num_cells
      call create_cell_locator(i, loc_p)
      call count_neighbours(loc_p, nnb)
      do j = 1, nnb
        call create_neighbour_locator(loc_p, j, loc_nb)
        call get_local_status(loc_nb, cell_local)
        if (cell_local) then
          call get_local_index(loc_nb, idx)
          if (ctr > size(adjncy)) then
            allocate (tmp(size(adjncy) + alloc_step))
            tmp(1:size(adjncy)) = adjncy(:)
            adjncy = tmp
            deallocate (tmp)
          end if
          adjncy(ctr) = idx
          ctr = ctr + 1
        end if
      end do
      xadj(i + 1) = ctr
    end do

    if (size(adjncy) > (ctr - 1)) then
      ! Shrink adjncy array
      allocate (tmp(ctr - 1))
      tmp(:) = adjncy(1:(ctr - 1))
      adjncy = tmp
    end if

  end subroutine build_adjacency_matrix

  !v Sets the offsets used for indexing into shared arrays for data that belongs to each rank.
  !  The halo cells may be interleaved with the local cells for some data so we need to store offsets
  !  for both types of arrays.
  subroutine set_offsets(shared_env, mesh)
    class(parallel_environment), intent(in) :: shared_env   !< The shared environment
    type(ccs_mesh), intent(inout) :: mesh                   !< The mesh

    integer(ccs_int), dimension(:), pointer :: shared_array_local_offsets   !< Offset within shared arrays for quantities that are locally indexed (i.e. each rank is responsible for local_num_cells of these)
    integer(ccs_int), dimension(:), pointer :: shared_array_total_offsets   !< Offset within shared arrays for quantities that are totally indexed (i.e. each rank is responsible for total_num_cells of these)
    integer :: shared_array_local_offsets_window                             !< Assoicated shared window
    integer :: shared_array_total_offsets_window                             !< Assoicated shared window
    integer(ccs_int), dimension(:), allocatable :: temp_offset
    integer(ccs_int) :: rank
    integer(ccs_int) :: i

    select type (shared_env)
    type is (parallel_environment_mpi)
      call create_shared_array(shared_env, shared_env%num_procs, shared_array_local_offsets, shared_array_local_offsets_window)
      call create_shared_array(shared_env, shared_env%num_procs, shared_array_total_offsets, shared_array_total_offsets_window)

      rank = shared_env%proc_id

      shared_array_local_offsets(rank + 1) = mesh%topo%local_num_cells
      shared_array_total_offsets(rank + 1) = mesh%topo%total_num_cells

      call sync(shared_env)
      if (rank == 0) then
        allocate (temp_offset(shared_env%num_procs))
        temp_offset(1) = 0
        do i = 2, shared_env%num_procs
          temp_offset(i) = temp_offset(i - 1) + shared_array_local_offsets(i - 1)
        end do
        shared_array_local_offsets = temp_offset

        do i = 2, shared_env%num_procs
          temp_offset(i) = temp_offset(i - 1) + shared_array_total_offsets(i - 1)
        end do
        shared_array_total_offsets = temp_offset
      end if
      call sync(shared_env)

      mesh%topo%shared_array_local_offset = shared_array_local_offsets(rank + 1)
      mesh%topo%shared_array_total_offset = shared_array_total_offsets(rank + 1)

      call destroy_shared_array(shared_env, shared_array_local_offsets, shared_array_local_offsets_window)
      call destroy_shared_array(shared_env, shared_array_total_offsets, shared_array_total_offsets_window)
    class default
      call error_abort("invalid parallel environment")
    end select
  end subroutine set_offsets

  !v check if parent of loc_nb (neighbour locator) is in loc_nb neighbour list, fail using assert if not
  subroutine test_mesh_internal_neighbours(loc_nb)

    use testing_lib, only: assert_neq, assert_bool

    type(neighbour_locator), intent(in) :: loc_nb ! Neighbour locator to test

    integer(ccs_int) :: index_nb
    integer(ccs_int) :: nnb
    integer(ccs_int) :: j
    logical :: found_parent
    logical :: is_boundary
    logical :: is_local
    type(cell_locator) :: cell_loc_nb
    type(neighbour_locator) :: loc_nb_nb

    associate (parent_idx => loc_nb%index_p)
      call get_local_index(loc_nb, index_nb)

      ! Neighbour index should not be its parents
      call assert_neq(index_nb, parent_idx, "Neighbour has same index as parent cell")

      call get_local_status(loc_nb, is_local)
      if (is_local) then
        ! Parent should be in neighbour's neighbour list
        call create_cell_locator(index_nb, cell_loc_nb)
        call count_neighbours(cell_loc_nb, nnb)
        found_parent = .false.
        do j = 1, nnb
          call create_neighbour_locator(cell_loc_nb, j, loc_nb_nb)
          call get_boundary_status(loc_nb_nb, is_boundary)
          if (.not. is_boundary) then ! We are looking for parent cell - by definition not a boundary!
            call get_local_index(loc_nb_nb, index_nb)
            if (index_nb == parent_idx) then
              found_parent = .true.
              exit
            end if
          end if
        end do
        call assert_bool(found_parent, "Couldn't find cell in neighbour's neighbour list")
      end if

    end associate

  end subroutine test_mesh_internal_neighbours

  subroutine report_mesh_surface_integral(par_env, mesh)
    use kinds, only: CCS_MPI_PRECISION

    class(parallel_environment), intent(in) :: par_env
    type(ccs_mesh), intent(inout) :: mesh

    type(face_locator) :: loc_f
    real(ccs_real), dimension(ndim) :: normal
    real(ccs_real) :: cell_integral
    real(ccs_real) :: area
    real(ccs_real) :: total_local, total_global
    real(ccs_real) :: avg_global
    real(ccs_real) :: max_local, max_global
    integer(ccs_int) :: n_cells
    integer(ccs_int) :: i, j
    integer :: ierr

    total_local = 0.0_ccs_real
    max_local = 0.0_ccs_real

    call set_mesh_object(mesh)
    call get_local_num_cells(n_cells)

    ! Compute the surface integral for each cell owned by this rank.
    do i = 1, n_cells
      cell_integral = 0.0_ccs_real

      do j = 1, mesh%topo%num_nb(i)
        call create_face_locator(i, j, loc_f)
        call get_face_area(loc_f, area)
        call get_face_normal(loc_f, normal)

        cell_integral = cell_integral + sum(normal(:)) * area
      end do

      ! Accumulate stats
      total_local = total_local + cell_integral
      max_local = max(max_local, cell_integral)

    end do

    avg_global = 0.0_ccs_real

    call nullify_mesh_object()

    ! Reductions to get global mesh diagnostics
    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_Allreduce(total_local, total_global, 1, CCS_MPI_PRECISION, MPI_SUM, par_env%comm, ierr)
      call MPI_Allreduce(max_local, max_global, 1, CCS_MPI_PRECISION, MPI_MAX, par_env%comm, ierr)
    class default
      call error_abort("invalid parallel environment")
    end select

    avg_global = total_global / real(mesh%topo%global_num_cells, ccs_real)

    if (is_root(par_env)) then
      write(log_unit_out,*) "* Total cell surface integral: ", total_global
      write(log_unit_out,*) "* Average cell surface integral: ", avg_global
      write(log_unit_out,*) "* Max cell surface integral: ", max_global
      write(log_unit_out,*) "* Max cell surface integral / average cell surface integral: ", max_global / avg_global
    end if

  end subroutine report_mesh_surface_integral

end module mesh_utils
