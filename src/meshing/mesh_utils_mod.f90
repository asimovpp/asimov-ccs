module mesh_utils
#include "ccs_macros.inc"

  use constants, only: ndim, geoext, adiosconfig
  use utils, only: exit_print, str
  use kinds, only: ccs_int, ccs_long, ccs_real
  use types, only: ccs_mesh, topology, geometry, &
                   io_environment, io_process, &
                   face_locator, cell_locator, neighbour_locator, vert_locator
  use io, only: read_scalar, read_array, &
                write_scalar, write_array, &
                configure_io, open_file, close_file, &
                initialise_io, cleanup_io
  use parallel, only: read_command_line_arguments
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi
  use meshing, only: get_global_index, get_local_index, count_neighbours, &
                     set_cell_location, set_neighbour_location, set_face_location, set_vert_location, &
                     set_face_index, get_boundary_status, get_local_status, &
                     get_centre, set_centre, &
                     set_area, set_normal
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

  integer, parameter :: front_bottom_left = 1_ccs_int
  integer, parameter :: front_bottom_right = 2_ccs_int
  integer, parameter :: front_top_right = 3_ccs_int
  integer, parameter :: front_top_left = 4_ccs_int
  integer, parameter :: back_bottom_left = 5_ccs_int
  integer, parameter :: back_bottom_right = 6_ccs_int
  integer, parameter :: back_top_right = 7_ccs_int
  integer, parameter :: back_top_left = 8_ccs_int

  private
  public :: build_square_mesh
  public :: build_mesh
  public :: read_mesh
  public :: write_mesh
  public :: global_start
  public :: local_count
  public :: count_mesh_faces
  public :: set_cell_face_indices
  public :: print_topo
  public :: print_geo

contains

  !v Read mesh from file
  subroutine read_mesh(par_env, case_name, mesh)

  use partitioning, only: partition_kway, compute_connectivity, compute_partitioner_input

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

    geo_file = case_name // "_mesh"// geoext
    adios2_file = case_name // adiosconfig

    call initialise_io(par_env, adios2_file, io_env)
    call configure_io(io_env, "geo_reader", geo_reader)

    call open_file(geo_file, "read", geo_reader)

    call read_topology(par_env, geo_reader, mesh)

    call compute_partitioner_input(par_env, mesh)
    call partition_kway(par_env, mesh)

    !call partition_kway(par_env, mesh)
    call partition_stride(par_env, mesh)

    ! for debug: Hardcode a 2 rank partitioning (even in serial)
    !mesh%topo%global_partition(1:mesh%topo%global_num_cells/2) = 1_ccs_int
    !mesh%topo%global_partition(mesh%topo%global_num_cells/2+1:mesh%topo%global_num_cells) = 0_ccs_int
    
    call compute_connectivity(par_env, mesh)

    call read_geometry(geo_reader, mesh)

    !call print_topo(par_env, mesh)
    !call print_geo(par_env, mesh)

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
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    class(io_process) :: geo_reader                                         !< The IO process for reading the file
    type(ccs_mesh), intent(inout) :: mesh                                   !< The mesh that will be read

    integer(ccs_int) :: i, j, k
    integer(ccs_int) :: num_bnd !< global number of boundary faces
    integer(ccs_int), dimension(:), allocatable :: bnd_rid, bnd_face
    integer(ccs_long), dimension(1) :: sel_start
    integer(ccs_long), dimension(1) :: sel_count

    integer(ccs_long), dimension(2) :: sel2_start
    integer(ccs_long), dimension(2) :: sel2_count

    ! Read attribute "ncel" - the total number of cells
    call read_scalar(geo_reader, "ncel", mesh%topo%global_num_cells)
    ! Read attribute "nfac" - the total number of faces
    call read_scalar(geo_reader, "nfac", mesh%topo%global_num_faces)
    ! Read attribute "maxfaces" - the maximum number of faces per cell
    call read_scalar(geo_reader, "maxfaces", mesh%topo%max_faces)
    ! Read attribute "nvrt" - the total number of vertices
    call read_scalar(geo_reader, "nvrt", mesh%topo%global_num_vertices)

    ! Read attribute "nbnd" - the total number of boundary faces
    call read_scalar(geo_reader, "nbnd", num_bnd)

    if (mesh%topo%max_faces == 6) then ! if cell are hexes
      mesh%topo%vert_per_cell = 8 ! 8 vertices per cell
    else
      call error_abort("Currently only supporting hex cells.")
    end if

    allocate (mesh%topo%face_cell1(mesh%topo%global_num_faces))
    allocate (mesh%topo%face_cell2(mesh%topo%global_num_faces))
    allocate (mesh%topo%bnd_rid(mesh%topo%global_num_faces))
    allocate (bnd_rid(num_bnd))
    allocate (bnd_face(num_bnd))
    allocate (mesh%topo%global_face_indices(mesh%topo%max_faces, mesh%topo%global_num_cells))
    allocate (mesh%topo%global_vertex_indices(mesh%topo%vert_per_cell, mesh%topo%global_num_cells))

    sel_start(1) = 0 ! Global index to start reading from
    sel_count(1) = mesh%topo%global_num_faces ! How many elements to read in total

    ! Read arrays face/cell1 and face/cell2
    call read_array(geo_reader, "/face/cell1", sel_start, sel_count, mesh%topo%face_cell1)
    call read_array(geo_reader, "/face/cell2", sel_start, sel_count, mesh%topo%face_cell2)

    sel_start(1) = 0 ! Global index to start reading from
    sel_count(1) = num_bnd ! How many elements to read in total
    call read_array(geo_reader, "/bnd/rid", sel_start, sel_count, bnd_rid)
    call read_array(geo_reader, "/bnd/face", sel_start, sel_count, bnd_face)

    ! make sure inside faces=0 and boundary faces are negative
    mesh%topo%bnd_rid(:) = 0_ccs_int
    mesh%topo%bnd_rid(bnd_face(:)) = - (bnd_rid(:) + 1_ccs_int)
    
    sel2_start = 0
    sel2_count(1) = mesh%topo%max_faces! topo%global_num_cells
    sel2_count(2) = mesh%topo%global_num_cells

    call read_array(geo_reader, "/cell/cface", sel2_start, sel2_count, mesh%topo%global_face_indices)

    sel2_count(1) = mesh%topo%vert_per_cell

    call read_array(geo_reader, "/cell/vertices", sel2_start, sel2_count, mesh%topo%global_vertex_indices)

    ! Create and populate the vtxdist array based on the total number of cells
    ! and the total number of ranks in the parallel environment
    allocate (mesh%topo%vtxdist(par_env%num_procs + 1)) ! vtxdist array is of size num_procs + 1 on all ranks

    mesh%topo%vtxdist(1) = 1                                                  ! First element is 1
    mesh%topo%vtxdist(par_env%num_procs + 1) = mesh%topo%global_num_cells + 1 ! Last element is total number of cells + 1

    ! Divide the total number of cells by the world size to
    ! compute the chunk sizes
    k = int(real(mesh%topo%global_num_cells) / par_env%num_procs)
    j = 1

    do i = 1, par_env%num_procs
      mesh%topo%vtxdist(i) = j
      j = j + k
    end do


  end subroutine read_topology

  !v Read the geometry data from an input (HDF5) file
  subroutine read_geometry(geo_reader, mesh)

    class(io_process) :: geo_reader                                         !< The IO process for reading the file
    type(ccs_mesh), intent(inout) :: mesh                                   !< The mesh%geometry that will be read

    integer(ccs_int) :: i, j, n, global_icell, local_icell
    integer(ccs_int) :: ier
    integer(ccs_int) :: vert_per_cell

    integer(ccs_long), dimension(1) :: vol_p_start
    integer(ccs_long), dimension(1) :: vol_p_count
    integer(ccs_long), dimension(2) :: x_p_start
    integer(ccs_long), dimension(2) :: x_p_count
    integer(ccs_long), dimension(2) :: f_xn_start
    integer(ccs_long), dimension(2) :: f_xn_count
    integer(ccs_long), dimension(1) :: f_a_start
    integer(ccs_long), dimension(1) :: f_a_count

    real(ccs_real), dimension(:), allocatable :: temp_vol_c ! Temp array for cell volumes
    real(ccs_real), dimension(:, :), allocatable :: temp_x_p ! Temp array for face centres
    real(ccs_real), dimension(:, :), allocatable :: temp_x_f ! Temp array for face centres
    real(ccs_real), dimension(:, :), allocatable :: temp_n_f ! Temp array for face normals
    real(ccs_real), dimension(:, :), allocatable :: temp_v_c ! Temp array for vertex coordinates
    real(ccs_real), dimension(:), allocatable :: temp_a_f ! Temp array for face areas

    type(face_locator) :: loc_f ! Face locator object
    type(vert_locator) :: loc_v ! Vertex locator object
    
    if (mesh%topo%max_faces == 6) then ! if cell are hexes
      vert_per_cell = 8 ! 8 vertices per cell
    else
      call error_abort("Currently only supporting hex cells.")
    end if

    ! Read attribute "scalefactor"
    call read_scalar(geo_reader, "scalefactor", mesh%geo%scalefactor)

    ! Starting point for reading chunk of data
    vol_p_start = 0 
    ! How many data points will be read?
    vol_p_count = mesh%topo%global_num_cells

    ! Allocate memory for cell volumes array on each MPI rank
    allocate (mesh%geo%volumes(mesh%topo%total_num_cells))
    allocate (temp_vol_c(mesh%topo%global_num_cells))

    ! Read variable "/cell/vol"
    call read_array(geo_reader, "/cell/vol", vol_p_start, vol_p_count, temp_vol_c)
    mesh%geo%volumes(:) = temp_vol_c(mesh%topo%global_indices(:))

    ! Starting point for reading chunk of data
    x_p_start = (/0, 0/)
    ! How many data points will be read?
    x_p_count = (/ndim, mesh%topo%global_num_cells/)

    ! Allocate memory for cell centre coordinates array on each MPI rank
    allocate (mesh%geo%x_p(ndim, mesh%topo%total_num_cells))
    allocate (temp_x_p(ndim, mesh%topo%global_num_cells))

    ! Read variable "/cell/x"
    call read_array(geo_reader, "/cell/x", x_p_start, x_p_count, temp_x_p)
    mesh%geo%x_p(:, :) = temp_x_p(:, mesh%topo%global_indices(:)) 

    ! Allocate temporary arrays for face centres, face normals, face areas and vertex coords
    allocate (temp_x_f(ndim, mesh%topo%global_num_faces))
    allocate (temp_n_f(ndim, mesh%topo%global_num_faces))
    allocate (temp_v_c(ndim, mesh%topo%global_num_vertices))
    allocate (temp_a_f(mesh%topo%global_num_faces))

    f_xn_start = 0
    f_xn_count(1) = ndim
    f_xn_count(2) = mesh%topo%global_num_faces

    ! Read variable "/face/x"
    call read_array(geo_reader, "/face/x", f_xn_start, f_xn_count, temp_x_f)
    ! Read variable "/face/n"
    call read_array(geo_reader, "/face/n", f_xn_start, f_xn_count, temp_n_f)

    f_xn_count(1) = ndim
    f_xn_count(2) = mesh%topo%global_num_vertices

    ! Read variable "/vert"
    call read_array(geo_reader, "/vert", f_xn_start, f_xn_count, temp_v_c)

    f_a_start = 0
    f_a_count(1) = mesh%topo%global_num_faces

    ! Read variable "/face/area"
    call read_array(geo_reader, "/face/area", f_a_start, f_a_count, temp_a_f)

    ! Allocate arrays for face centres, face normals, face areas arrand vertex coordinates
    allocate (mesh%geo%x_f(ndim, mesh%topo%max_faces, mesh%topo%local_num_cells))
    allocate (mesh%geo%face_normals(ndim, mesh%topo%max_faces, mesh%topo%local_num_cells))
    allocate (mesh%geo%face_areas(mesh%topo%max_faces, mesh%topo%local_num_cells))
    allocate (mesh%geo%vert_coords(ndim, vert_per_cell, mesh%topo%local_num_cells))


    !    do k = start, end ! loop over cells owned by current process
    do local_icell = 1, mesh%topo%local_num_cells ! loop over cells owned by current process
      global_icell = mesh%topo%global_indices(local_icell)

      do j = 1, mesh%topo%max_faces ! loop over all faces for each cell
        call set_face_location(mesh, local_icell, j, loc_f)

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
        call set_vert_location(mesh, local_icell, j, loc_v)
         
        n = mesh%topo%global_vertex_indices(j, global_icell)
        call set_centre(loc_v, temp_v_c(:, n))
      end do

    end do

    ! Delete temp arrays
    deallocate (temp_vol_c, stat=ier)
    deallocate (temp_x_f, stat=ier)
    deallocate (temp_v_c, stat=ier)
    deallocate (temp_x_p, stat=ier)
    deallocate (temp_n_f, stat=ier)
    deallocate (temp_a_f, stat=ier)

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

    ! Set ADIOS2 config file name
    adios2_file = case_name // adiosconfig

    ! Set geo file name
    geo_file = case_name // geoext 

    ! Open geo file for writing
    call initialise_io(par_env, adios2_file, io_env)
    call configure_io(io_env, "geo_writer", geo_writer)
    call open_file(geo_file, "write", geo_writer)

    ! Write mesh
    call write_topology(geo_writer, mesh)
    call write_geometry(par_env, geo_writer, mesh)

    ! Close the file and ADIOS2 engine
    call close_file(geo_writer)

    ! Finalise the ADIOS2 environment
    call cleanup_io(io_env)

  end subroutine write_mesh

  !v Write the mesh topology data to file
  subroutine write_topology(geo_writer, mesh)

    ! Arguments
    class(io_process), allocatable, target, intent(in) :: geo_writer        !< The IO process for writing the mesh ("geo") file
    type(ccs_mesh), intent(in) :: mesh                                      !< The mesh

    ! Local variables
    integer(ccs_long), dimension(2) :: sel2_shape
    integer(ccs_long), dimension(2) :: sel2_start
    integer(ccs_long), dimension(2) :: sel2_count

    ! Write cell vertices
    sel2_shape(1) = mesh%topo%vert_per_cell
    sel2_shape(2) = mesh%topo%global_num_cells
    sel2_start(1) = 0
    sel2_start(2) = mesh%topo%global_indices(1) - 1
    sel2_count(1) = mesh%topo%vert_per_cell
    sel2_count(2) = mesh%topo%local_num_cells

    call write_array(geo_writer, "/cell/vertices", sel2_shape, sel2_start, sel2_count, mesh%topo%global_vertex_indices)

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

    real(ccs_real), dimension(:, :), allocatable :: vert_coords_tmp

    ! Root process calculates all (global) vertex coords and stores temporarily
    if (par_env%proc_id == par_env%root) then
      allocate (vert_coords_tmp(ndim, mesh%topo%global_num_vertices))

      vert_coords_tmp = 0.0_ccs_real

      if (mesh%topo%vert_per_cell == 4) then
        verts_per_side = nint(mesh%topo%global_num_vertices**(1./2))
      else
        verts_per_side = nint(mesh%topo%global_num_vertices**(1./3))
      end if

      do i = 1, mesh%topo%global_num_vertices
        vert_coords_tmp(1, i) = modulo(i - 1, verts_per_side) * mesh%geo%h
        vert_coords_tmp(2, i) = modulo((i - 1) / verts_per_side, verts_per_side) * mesh%geo%h
        vert_coords_tmp(3, i) = ((i - 1) / (verts_per_side * verts_per_side)) * mesh%geo%h
      end do
    end if

    sel2_shape(1) = ndim
    sel2_shape(2) = mesh%topo%global_num_vertices
    sel2_start(1) = 0
    sel2_start(2) = 0
    if (par_env%proc_id == par_env%root) then
      sel2_count(1) = ndim
      sel2_count(2) = mesh%topo%global_num_vertices
    else
      sel2_count(1) = 0
      sel2_count(2) = 0
    end if

    call write_array(geo_writer, "/vert", sel2_shape, sel2_start, sel2_count, vert_coords_tmp)

    if (par_env%proc_id == par_env%root) then
      deallocate (vert_coords_tmp)
    end if

  end subroutine write_geometry

 !v Utility constructor to build a 3D mesh with hex cells.
  !
  !  Builds a Cartesian grid of nx*ny*nz cells.
  function build_square_mesh(par_env, cps, side_length) result(mesh)

  use partitioning, only: partition_kway, compute_connectivity, compute_partitioner_input

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    integer(ccs_int), intent(in) :: cps                !< Number of cells per side of the mesh.
    real(ccs_real), intent(in) :: side_length          !< The length of the side.

    type(ccs_mesh) :: mesh                             !< The resulting mesh.

    call build_square_topology(par_env, cps, mesh)

    call compute_partitioner_input(par_env, mesh)

    !call partition_kway(par_env, mesh)
    call partition_stride(par_env, mesh)

    ! for debug: Hardcode a 2 rank partitioning (even in serial)
    !mesh%topo%global_partition(1:mesh%topo%global_num_cells/2) = 1_ccs_int
    !mesh%topo%global_partition(mesh%topo%global_num_cells/2+1:mesh%topo%global_num_cells) = 0_ccs_int
    
    call compute_connectivity(par_env, mesh)

    call build_square_geometry(par_env, cps, side_length, mesh)

  end function build_square_mesh


  !v Utility constructor to build a square mesh.
  !
  !  Builds a Cartesian grid of NxN cells on the domain LxL.
  subroutine build_square_topology(par_env, cps, mesh)

    class(parallel_environment), intent(in) :: par_env !< The parallel environment to construct the mesh.
    integer(ccs_int), intent(in) :: cps                !< Number of cells per side of the mesh.

    type(ccs_mesh), intent(inout) :: mesh                             !< The resulting mesh.

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
      mesh%topo%global_num_vertices = (cps + 1)**2

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

        ! Set number of vertices per cell
        mesh%topo%vert_per_cell = 4_ccs_int

        ! Allocate mesh arrays
        allocate (mesh%topo%global_indices(mesh%topo%local_num_cells))
        allocate (mesh%topo%num_nb(mesh%topo%local_num_cells))
        allocate (mesh%topo%nb_indices(mesh%topo%max_faces, mesh%topo%local_num_cells))
        allocate (mesh%topo%face_indices(mesh%topo%max_faces, mesh%topo%local_num_cells))
        allocate (mesh%topo%global_vertex_indices(mesh%topo%vert_per_cell, mesh%topo%local_num_cells))

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

      ! Set the global mesh parameters
      mesh%topo%global_num_cells = cps**2
      mesh%topo%global_num_vertices = (cps + 1)**2
      mesh%topo%num_faces = count_mesh_faces(mesh)

      call set_cell_face_indices(mesh)

    class default
      call error_abort("Unknown parallel environment type.")

    end select

   end subroutine build_square_topology



  subroutine build_square_geometry(par_env, cps, side_length, mesh)
  
    class(parallel_environment), intent(in) :: par_env !< The parallel environment to construct the mesh.
    integer(ccs_int), intent(in) :: cps                !< Number of cells per side of the mesh.
    real(ccs_real), intent(in) :: side_length          !< The length of each side.

    type(ccs_mesh), intent(inout) :: mesh                             !< The resulting mesh.

    integer(ccs_int) :: i               ! Loop counter
    integer(ccs_int) :: ii              ! Zero-indexed loop counter (simplifies some operations)
    integer(ccs_int) :: face_counter    ! Cell-local face counter
    integer(ccs_int) :: vertex_counter  ! Cell-local vertex counter

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

      if (allocated(mesh%topo%global_vertex_indices)) then
        deallocate (mesh%topo%global_vertex_indices)
      end if
      allocate (mesh%topo%global_vertex_indices(mesh%topo%vert_per_cell, mesh%topo%local_num_cells))

      ! Global vertex numbering
      do i = 1, mesh%topo%local_num_cells
        associate (global_vert_index => mesh%topo%global_vertex_indices(:, i))
          ii = mesh%topo%global_indices(i)

          global_vert_index(front_bottom_left) = ii + (ii - 1) / cps
          global_vert_index(front_bottom_right) = ii + (ii - 1) / cps + 1
          global_vert_index(front_top_left) = ii + cps + (ii + cps - 1) / cps
          global_vert_index(front_top_right) = ii + cps + (ii + cps - 1) / cps + 1
        end associate
      end do


      allocate (mesh%geo%x_p(ndim, mesh%topo%total_num_cells))
      allocate (mesh%geo%x_f(ndim, mesh%topo%max_faces, mesh%topo%local_num_cells)) !< @note Currently hardcoded as a 2D mesh. @endnote
      allocate (mesh%geo%volumes(mesh%topo%total_num_cells))
      allocate (mesh%geo%face_areas(mesh%topo%max_faces, mesh%topo%local_num_cells))
      allocate (mesh%geo%face_normals(ndim, mesh%topo%max_faces, mesh%topo%local_num_cells)) ! Currently hardcoded as a 2D mesh.
      allocate (mesh%geo%vert_coords(ndim, mesh%topo%vert_per_cell, mesh%topo%local_num_cells))

      mesh%geo%h = side_length / real(cps, ccs_real)
      mesh%geo%volumes(:) = mesh%geo%h**2 !< @note Mesh is square and 2D @endnote
      mesh%geo%face_normals(:, :, :) = 0.0_ccs_real
      mesh%geo%x_p(:, :) = 0.0_ccs_real
      mesh%geo%x_f(:, :, :) = 0.0_ccs_real
      mesh%geo%face_areas(:, :) = mesh%geo%h  ! Mesh is square and 2D
      mesh%geo%vert_coords(:, :, :) = 0.0_ccs_real

      ! Set cell centre
      associate (h => mesh%geo%h)
        do i = 1_ccs_int, mesh%topo%total_num_cells
          call set_cell_location(mesh, i, loc_p)
           
          ii = mesh%topo%global_indices(i)

          x_p(1) = (modulo(ii - 1, cps) + 0.5_ccs_real) * h
          x_p(2) = ((ii - 1) / cps + 0.5_ccs_real) * h

          call set_centre(loc_p, x_p)
        end do

        do i = 1_ccs_int, mesh%topo%local_num_cells
          call set_cell_location(mesh, i, loc_p)
          call get_centre(loc_p, x_p)
          

          do face_counter = 1_ccs_int, mesh%topo%max_faces
            
            call set_neighbour_location(loc_p, face_counter, loc_nb)
            call get_boundary_status(loc_nb, is_boundary)
            call set_face_location(mesh, i, face_counter, loc_f)

            if (.not. is_boundary) then
              ! faces are midway between cell centre and nb cell centre
              call get_centre(loc_nb, x_nb_3)
              x_nb(:) = x_nb_3(1:2) ! XXX: hacky fix for issue with resolving get_neighbour_centre in 2D

              x_f(:) = 0.5_ccs_real*(x_p(:) + x_nb(:))
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

        do i = 1_ccs_int, mesh%topo%local_num_cells
          call set_cell_location(mesh, i, loc_p)
          call get_centre(loc_p, x_p)

          vertex_counter = front_bottom_left
          call set_vert_location(mesh, i, vertex_counter, loc_v)
          x_v(1) = x_p(1) - 0.5_ccs_real * h
          x_v(2) = x_p(2) - 0.5_ccs_real * h
          call set_centre(loc_v, x_v)
          
          vertex_counter = front_bottom_right
          call set_vert_location(mesh, i, vertex_counter, loc_v)
          x_v(1) = x_p(1) + 0.5_ccs_real * h
          x_v(2) = x_p(2) - 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = front_top_left
          call set_vert_location(mesh, i, vertex_counter, loc_v)
          x_v(1) = x_p(1) - 0.5_ccs_real * h
          x_v(2) = x_p(2) + 0.5_ccs_real * h
          call set_centre(loc_v, x_v)

          vertex_counter = front_top_right
          call set_vert_location(mesh, i, vertex_counter, loc_v)
          x_v(1) = x_p(1) + 0.5_ccs_real * h
          x_v(2) = x_p(2) + 0.5_ccs_real * h
          call set_centre(loc_v, x_v)
        end do
      end associate

    class default
      call error_abort("Unknown parallel environment type.")

    end select

  end subroutine build_square_geometry

  !v Utility constructor to build a 3D mesh with hex cells.
  !
  !  Builds a Cartesian grid of nx*ny*nz cells.
  function build_mesh(par_env, nx, ny, nz, side_length) result(mesh)

  use partitioning, only: partition_kway, compute_connectivity, compute_partitioner_input

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    integer(ccs_int), intent(in) :: nx                 !< Number of cells in the x direction.
    integer(ccs_int), intent(in) :: ny                 !< Number of cells in the y direction.
    integer(ccs_int), intent(in) :: nz                 !< Number of cells in the z direction.
    real(ccs_real), intent(in) :: side_length          !< The length of the side.

    type(ccs_mesh) :: mesh                             !< The resulting mesh.

    if (.not. (nx .eq. ny .and. ny .eq. nz)) then !< @note Must be a cube (for now) @endnote
      call error_abort("Only supporting cubes for now - nx, ny and nz must be the same!")
    end if

    call build_topology(par_env, nx, ny, nz, mesh)

    call compute_partitioner_input(par_env, mesh)

    !call partition_kway(par_env, mesh)
    call partition_stride(par_env, mesh)

    ! for debug: Hardcode a 2 rank partitioning (even in serial)
    !mesh%topo%global_partition(1:mesh%topo%global_num_cells/2) = 1_ccs_int
    !mesh%topo%global_partition(mesh%topo%global_num_cells/2+1:mesh%topo%global_num_cells) = 0_ccs_int
    
    call compute_connectivity(par_env, mesh)

    call build_geometry(par_env, nx, ny, nz, side_length, mesh)

  end function build_mesh
 
  !v Utility constructor to build a 3D mesh with hex cells.
  !
  !  Builds a Cartesian grid of nx*ny*nz cells.
  subroutine build_topology(par_env, nx, ny, nz, mesh)

    class(parallel_environment), intent(in) :: par_env !< The parallel environment to construct the mesh.
    integer(ccs_int), intent(in) :: nx                 !< Number of cells in the x direction.
    integer(ccs_int), intent(in) :: ny                 !< Number of cells in the y direction.
    integer(ccs_int), intent(in) :: nz                 !< Number of cells in the z direction.

    type(ccs_mesh), intent(inout) :: mesh                             !< The resulting mesh.

    integer(ccs_int) :: start_global    ! The (global) starting index of a partition
    integer(ccs_int) :: end_global      ! The (global) last index of a partition
    integer(ccs_int) :: i,j,k           ! Loop counter
    integer(ccs_int) :: ii              ! Zero-indexed loop counter (simplifies some operations)
    integer(ccs_int) :: index_counter   ! Local index counter
    integer(ccs_int) :: face_counter    ! Cell-local face counter
    integer(ccs_int) :: face_index_counter    ! global face counter

    integer(ccs_int) :: index_nb        ! The local index of a neighbour cell
    integer(ccs_int) :: global_index_nb ! The global index of a neighbour cell

    select type (par_env)
    type is (parallel_environment_mpi)

      ! Set the global mesh parameters
      mesh%topo%global_num_cells = nx * ny * nz
      mesh%topo%global_num_vertices = (nx + 1) * (ny + 1) * (nz + 1)

      ! Associate aliases to make code easier to read
      associate (nglobal => mesh%topo%global_num_cells)

        ! Determine ownership range
        start_global = global_start(nglobal, par_env%proc_id, par_env%num_procs)
        mesh%topo%local_num_cells = local_count(nglobal, par_env%proc_id, par_env%num_procs)
        end_global = start_global + (mesh%topo%local_num_cells - 1)

        ! Set max number of faces (constant, 6)
        mesh%topo%max_faces = 6_ccs_int

        ! Set number of vertices per cell (constant, 8)
        mesh%topo%vert_per_cell = 8

        ! Allocate mesh arrays
        allocate (mesh%topo%global_indices(mesh%topo%local_num_cells))
        allocate (mesh%topo%num_nb(mesh%topo%local_num_cells))
        allocate (mesh%topo%nb_indices(mesh%topo%max_faces, mesh%topo%local_num_cells))
        allocate (mesh%topo%face_indices(mesh%topo%max_faces, mesh%topo%local_num_cells))
        allocate (mesh%topo%global_vertex_indices(mesh%topo%vert_per_cell, mesh%topo%local_num_cells))

        ! Initialise mesh arrays
        mesh%topo%num_nb(:) = mesh%topo%max_faces ! All cells have 6 neighbours (possibly ghost/boundary cells)

        ! Initalise neighbour indices
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
          if (modulo(ii / nx, ny) == 0_ccs_int) then
            index_nb = -bottom
            global_index_nb = -bottom
          else
            index_nb = index_counter - nx
            global_index_nb = i - nx
          end if
          call build_local_mesh_add_neighbour(index_counter, face_counter, index_nb, global_index_nb, mesh)

          ! Construct top (4) face/neighbour
          face_counter = top
          if (modulo(ii / nx, ny) == (ny - 1_ccs_int)) then
            index_nb = -top
            global_index_nb = -top
          else
            index_nb = index_counter + nx
            global_index_nb = i + nx
          end if
          call build_local_mesh_add_neighbour(index_counter, face_counter, index_nb, global_index_nb, mesh)

          ! Construct back (5) face/neighbour
          face_counter = back
          if ((ii / (nx * ny)) == 0_ccs_int) then
            index_nb = -back
            global_index_nb = -back
          else
            index_nb = index_counter - nx * ny
            global_index_nb = i - nx * ny
          end if
          call build_local_mesh_add_neighbour(index_counter, face_counter, index_nb, global_index_nb, mesh)

          ! Construct front (6) face/neighbour
          face_counter = front
          if ((ii / (nx * ny)) == nz - 1_ccs_int) then
            index_nb = -front
            global_index_nb = -front
          else
            index_nb = index_counter + nx * ny
            global_index_nb = i + nx * ny
          end if
          call build_local_mesh_add_neighbour(index_counter, face_counter, index_nb, global_index_nb, mesh)

          index_counter = index_counter + 1_ccs_int

        end do
      end associate

      ! print*,"Neighbour indices: ",mesh%neighbour_indices

      mesh%topo%total_num_cells = size(mesh%topo%global_indices)
      mesh%topo%halo_num_cells = mesh%topo%total_num_cells - mesh%topo%local_num_cells
      

      mesh%topo%global_num_faces = (nx+1)*ny*nz + nx*(ny+1)*nz + nx*ny*(nz+1)
      allocate (mesh%topo%face_cell1(mesh%topo%global_num_faces))
      allocate (mesh%topo%face_cell2(mesh%topo%global_num_faces))
      allocate (mesh%topo%global_face_indices(mesh%topo%max_faces, mesh%topo%global_num_cells))
      allocate (mesh%topo%bnd_rid(mesh%topo%global_num_faces))

      mesh%topo%face_cell1(:) = 0_ccs_int
      mesh%topo%face_cell2(:) = 0_ccs_int
      mesh%topo%global_face_indices(:,:) = 0_ccs_int
      mesh%topo%bnd_rid(:) = 0_ccs_int

      ! Construct face_cell1 and face_cell2 following:
      !  - face_cell1 < face_cell2
      !  - and if face is a boundary, then: face_cell1 = current_cell, face_cell2 = 0
      face_index_counter = 1_ccs_int
      do i = 1, mesh%topo%global_num_cells

        ii = i - 1_ccs_int

        ! Construct left (1) face/neighbour
        face_counter = left
        if (modulo(ii, nx) == 0_ccs_int) then
          global_index_nb = -left
          mesh%topo%face_cell1(face_index_counter) = i
          mesh%topo%face_cell2(face_index_counter) = 0

          mesh%topo%global_face_indices(face_counter, i) = face_index_counter
          mesh%topo%bnd_rid(face_index_counter) = global_index_nb
          face_index_counter = face_index_counter + 1_ccs_int
        else
          !global_index_nb = i - 1_ccs_int
          !mesh%topo%face_cell1(face_index_counter) = global_index_nb
          !mesh%topo%face_cell2(face_index_counter) = i
        end if

        ! Construct right (2) face/neighbour
        face_counter = right
        mesh%topo%global_face_indices(face_counter, i) = face_index_counter
        if (modulo(ii, nx) == (nx - 1_ccs_int)) then
          global_index_nb = -right
          mesh%topo%face_cell1(face_index_counter) = i
          mesh%topo%face_cell2(face_index_counter) = 0 

          mesh%topo%global_face_indices(face_counter, i) = face_index_counter
          mesh%topo%bnd_rid(face_index_counter) = global_index_nb
        else
          global_index_nb = i + 1_ccs_int
          mesh%topo%face_cell1(face_index_counter) = i
          mesh%topo%face_cell2(face_index_counter) = global_index_nb

          mesh%topo%global_face_indices(face_counter, i) = face_index_counter
          mesh%topo%global_face_indices(left, global_index_nb) = face_index_counter
        end if
        face_index_counter = face_index_counter + 1_ccs_int

        ! Construct bottom (3) face/neighbour
        face_counter = bottom
        if (modulo(ii / nx, ny) == 0_ccs_int) then
          global_index_nb = -bottom
          mesh%topo%face_cell1(face_index_counter) = i
          mesh%topo%face_cell2(face_index_counter) = 0

          mesh%topo%global_face_indices(face_counter, i) = face_index_counter
          mesh%topo%bnd_rid(face_index_counter) = global_index_nb
          face_index_counter = face_index_counter + 1_ccs_int
        else
          !global_index_nb = i - nx
          !mesh%topo%face_cell1(face_index_counter) = global_index_nb
          !mesh%topo%face_cell2(face_index_counter) = i
        end if

        ! Construct top (4) face/neighbour
        face_counter = top
        if (modulo(ii / nx, ny) == (ny - 1_ccs_int)) then
          global_index_nb = -top
          mesh%topo%face_cell1(face_index_counter) = i
          mesh%topo%face_cell2(face_index_counter) = 0

          mesh%topo%global_face_indices(face_counter, i) = face_index_counter
          mesh%topo%bnd_rid(face_index_counter) = global_index_nb
        else
          global_index_nb = i + nx
          mesh%topo%face_cell1(face_index_counter) = i
          mesh%topo%face_cell2(face_index_counter) = global_index_nb

          mesh%topo%global_face_indices(face_counter, i) = face_index_counter
          mesh%topo%global_face_indices(bottom, global_index_nb) = face_index_counter
        end if
        face_index_counter = face_index_counter + 1_ccs_int

        ! Construct back (5) face/neighbour
        face_counter = back
        if ((ii / (nx * ny)) == 0_ccs_int) then
          global_index_nb = -back
          mesh%topo%face_cell1(face_index_counter) = i
          mesh%topo%face_cell2(face_index_counter) = 0

          mesh%topo%global_face_indices(face_counter, i) = face_index_counter
          mesh%topo%bnd_rid(face_index_counter) = global_index_nb
          face_index_counter = face_index_counter + 1_ccs_int
        else
          !global_index_nb = i - nx * ny
          !mesh%topo%face_cell1(face_index_counter) = global_index_nb
          !mesh%topo%face_cell2(face_index_counter) = i
        end if

        ! Construct front (6) face/neighbour
        face_counter = front
        if ((ii / (nx * ny)) == nz - 1_ccs_int) then
          global_index_nb = -front
          mesh%topo%face_cell1(face_index_counter) = i
          mesh%topo%face_cell2(face_index_counter) = 0

          mesh%topo%global_face_indices(face_counter, i) = face_index_counter
          mesh%topo%bnd_rid(face_index_counter) = global_index_nb
        else
          global_index_nb = i + nx * ny
          mesh%topo%face_cell1(face_index_counter) = i
          mesh%topo%face_cell2(face_index_counter) = global_index_nb

          mesh%topo%global_face_indices(face_counter, i) = face_index_counter
          mesh%topo%global_face_indices(back, global_index_nb) = face_index_counter
        end if
        face_index_counter = face_index_counter + 1_ccs_int

      end do
        

      mesh%topo%num_faces = count_mesh_faces(mesh)

      call set_cell_face_indices(mesh)


      ! Create and populate the vtxdist array based on the total number of cells
      ! and the total number of ranks in the parallel environment
      allocate (mesh%topo%vtxdist(par_env%num_procs + 1)) ! vtxdist array is of size num_procs + 1 on all ranks

      mesh%topo%vtxdist(1) = 1                                                  ! First element is 1
      mesh%topo%vtxdist(par_env%num_procs + 1) = mesh%topo%global_num_cells + 1 ! Last element is total number of cells + 1

      ! Divide the total number of cells by the world size to
      ! compute the chunk sizes
      k = int(real(mesh%topo%global_num_cells) / par_env%num_procs)
      j = 1

      do i = 1, par_env%num_procs
        mesh%topo%vtxdist(i) = j
        j = j + k
      end do
  
    class default
      call error_abort("Unknown parallel environment type.")

    end select

  end subroutine build_topology


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
    integer(ccs_int) :: a, b, c, d, e       ! Temporary variables

    logical :: is_boundary

    integer(ccs_int) :: index_nb        ! The local index of a neighbour cell

    real(ccs_real), dimension(3) :: x_p ! Cell centre array
    type(cell_locator) :: loc_p         ! Cell locator object

    real(ccs_real), dimension(3) :: x_nb ! Cell centre array of neighbour cell
    type(neighbour_locator) :: loc_nb    ! the neighbour locator object.

    real(ccs_real), dimension(3) :: x_f    ! Face centre array
    real(ccs_real), dimension(3) :: normal ! Face normal array
    type(face_locator) :: loc_f            ! Face locator object

    real(ccs_real), dimension(3) :: x_v ! Vertex centre array
    type(vert_locator) :: loc_v         ! Vertex locator object
 
    e = nz ! silence dummy argument unused warning

    select type (par_env)
    type is (parallel_environment_mpi)

      if (allocated(mesh%topo%global_vertex_indices)) then
        deallocate (mesh%topo%global_vertex_indices)
      end if
      allocate (mesh%topo%global_vertex_indices(mesh%topo%vert_per_cell, mesh%topo%local_num_cells))

      ! Global vertex numbering
      do i = 1, mesh%topo%local_num_cells
        associate (global_vert_index => mesh%topo%global_vertex_indices(:, i))
          ii = mesh%topo%global_indices(i)
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


      allocate (mesh%geo%x_p(ndim, mesh%topo%total_num_cells))
      allocate (mesh%geo%x_f(ndim, mesh%topo%max_faces, mesh%topo%local_num_cells))
      allocate (mesh%geo%volumes(mesh%topo%total_num_cells))
      allocate (mesh%geo%face_areas(mesh%topo%max_faces, mesh%topo%local_num_cells))
      allocate (mesh%geo%face_normals(ndim, mesh%topo%max_faces, mesh%topo%local_num_cells))
      allocate (mesh%geo%vert_coords(ndim, mesh%topo%vert_per_cell, mesh%topo%local_num_cells))

      mesh%geo%h = side_length / real(nx, ccs_real) !< @note Assumes cube @endnote
      mesh%geo%volumes(:) = mesh%geo%h**3 !< @note Mesh is cube @endnote
      mesh%geo%face_normals(:, :, :) = 0.0_ccs_real
      mesh%geo%x_p(:, :) = 0.0_ccs_real
      mesh%geo%x_f(:, :, :) = 0.0_ccs_real
      mesh%geo%face_areas(:, :) = mesh%geo%h**2
      mesh%geo%vert_coords(:, :, :) = 0.0_ccs_real

      ! Set cell centre
      associate (h => mesh%geo%h)
        do i = 1_ccs_int, mesh%topo%total_num_cells
          call set_cell_location(mesh, i, loc_p)
            
          ii = mesh%topo%global_indices(i)
          x_p(1) = (modulo(ii - 1, nx) + 0.5_ccs_real) * h
          x_p(2) = (modulo((ii - 1) / nx, ny) + 0.5_ccs_real) * h
          x_p(3) = (((ii - 1) / (nx * ny)) + 0.5_ccs_real) * h

          call set_centre(loc_p, x_p)
        end do

        do i = 1_ccs_int, mesh%topo%local_num_cells
          call set_cell_location(mesh, i, loc_p)
          call get_centre(loc_p, x_p)
          
          do face_counter = 1_ccs_int, mesh%topo%max_faces
            
            call set_neighbour_location(loc_p, face_counter, loc_nb)
            call get_boundary_status(loc_nb, is_boundary)
            call set_face_location(mesh, i, face_counter, loc_f)

            if (.not. is_boundary) then
              ! faces are midway between cell centre and nb cell centre
              call get_centre(loc_nb, x_nb)

              x_f(:) = 0.5_ccs_real*(x_p(:) + x_nb(:))
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

        do i = 1_ccs_int, mesh%topo%local_num_cells
            call set_cell_location(mesh, i, loc_p)
            call get_centre(loc_p, x_p)
          
            vertex_counter = front_bottom_left
            call set_vert_location(mesh, i, vertex_counter, loc_v)
            x_v(1) = x_p(1) - 0.5_ccs_real * h
            x_v(2) = x_p(2) - 0.5_ccs_real * h
            x_v(3) = x_p(3) + 0.5_ccs_real * h
            call set_centre(loc_v, x_v)
            
            vertex_counter = front_bottom_right
            call set_vert_location(mesh, i, vertex_counter, loc_v)
            x_v(1) = x_p(1) + 0.5_ccs_real * h
            x_v(2) = x_p(2) - 0.5_ccs_real * h
            x_v(3) = x_p(3) + 0.5_ccs_real * h
            call set_centre(loc_v, x_v)

            vertex_counter = front_top_left
            call set_vert_location(mesh, i, vertex_counter, loc_v)
            x_v(1) = x_p(1) - 0.5_ccs_real * h
            x_v(2) = x_p(2) + 0.5_ccs_real * h
            x_v(3) = x_p(3) + 0.5_ccs_real * h
            call set_centre(loc_v, x_v)

            vertex_counter = front_top_right
            call set_vert_location(mesh, i, vertex_counter, loc_v)
            x_v(1) = x_p(1) + 0.5_ccs_real * h
            x_v(2) = x_p(2) + 0.5_ccs_real * h
            x_v(3) = x_p(3) + 0.5_ccs_real * h
            call set_centre(loc_v, x_v)
            
            vertex_counter = back_bottom_left
            call set_vert_location(mesh, i, vertex_counter, loc_v)
            x_v(1) = x_p(1) - 0.5_ccs_real * h
            x_v(2) = x_p(2) - 0.5_ccs_real * h
            x_v(3) = x_p(3) - 0.5_ccs_real * h
            call set_centre(loc_v, x_v)

            vertex_counter = back_bottom_right
            call set_vert_location(mesh, i, vertex_counter, loc_v)
            x_v(1) = x_p(1) + 0.5_ccs_real * h
            x_v(2) = x_p(2) - 0.5_ccs_real * h
            x_v(3) = x_p(3) - 0.5_ccs_real * h
            call set_centre(loc_v, x_v)

            vertex_counter = back_top_left
            call set_vert_location(mesh, i, vertex_counter, loc_v)
            x_v(1) = x_p(1) - 0.5_ccs_real * h
            x_v(2) = x_p(2) + 0.5_ccs_real * h
            x_v(3) = x_p(3) - 0.5_ccs_real * h
            call set_centre(loc_v, x_v)

            vertex_counter = back_top_right
            call set_vert_location(mesh, i, vertex_counter, loc_v)
            x_v(1) = x_p(1) + 0.5_ccs_real * h
            x_v(2) = x_p(2) + 0.5_ccs_real * h
            x_v(3) = x_p(3) - 0.5_ccs_real * h
            call set_centre(loc_v, x_v)
        end do
      end associate

    class default
      call error_abort("Unknown parallel environment type.")

    end select

  end subroutine build_geometry 

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
    ! dragons. If this were intent(out) it
    ! would be deallocated on entry!
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
    integer(ccs_int) :: face_counter      ! Face index counter
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


  ! Populate mesh%topo%global_partition with a split of cells in stride using global_start and local_count
  subroutine partition_stride(par_env, mesh)
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(ccs_mesh), intent(inout) :: mesh                             !< The resulting mesh.

    integer(ccs_int) :: iproc, start, end

    do iproc = 0, par_env%num_procs -1
      start = global_start(mesh%topo%global_num_cells, iproc, par_env%num_procs)
      end = start + local_count(mesh%topo%global_num_cells, iproc, par_env%num_procs) -1
      mesh%topo%global_partition(start:end) = iproc
    end do

  end subroutine partition_stride

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

  ! Print mesh geometry object
  subroutine print_geo(par_env, mesh)
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(ccs_mesh), intent(in) :: mesh       !< the mesh
    integer(ccs_int) :: i          ! loop counters
    integer(ccs_int) :: nb_elem = 10
    !integer(ccs_int) :: file_unit


    print *, "############################# Print Geometry ########################################"

    print *, par_env%proc_id, "h                  : ", mesh%geo%h
    print *, par_env%proc_id, "scalefactor        : ", mesh%geo%scalefactor
    print *, ""

    if (allocated(mesh%geo%volumes))    print *, par_env%proc_id, "volumes         : ", mesh%geo%volumes(1:nb_elem)

    print *, ""
    if (allocated(mesh%geo%face_areas))   then
      do i=1, nb_elem
        print *, par_env%proc_id, "face_areas(1:"// str(nb_elem/2) // ", " // str(i) //")", mesh%geo%face_areas(1:nb_elem/2,i)
      end do
    end if

    print *, ""
    if (allocated(mesh%geo%x_p))   then
      do i=1, nb_elem
        print *, par_env%proc_id, "x_p(1:"// str(nb_elem/2) // ", " // str(i) //")", mesh%geo%x_p(1:nb_elem/2,i)
      end do
    end if

    print *, ""
    if (allocated(mesh%geo%x_f))   then
      do i=1, nb_elem
        print *, par_env%proc_id, "x_f(1:"// str(nb_elem/2) // ", " // str(i) //")", mesh%geo%x_f(2, 1:nb_elem/2,i)
      end do
    end if

    print *, ""
    if (allocated(mesh%geo%face_normals))   then
      do i=1, nb_elem
        print *, par_env%proc_id, "face_normals(1:"// str(nb_elem/2) // ", " // str(i) //")", mesh%geo%face_normals(2, 1:nb_elem/2,i)
      end do
    end if

    print *, ""
    if (allocated(mesh%geo%vert_coords))   then
      do i=1, nb_elem
        print *, par_env%proc_id, "vert_coords(1:"// str(nb_elem/2) // ", " // str(i) //")", mesh%geo%vert_coords(2, 1:nb_elem/2,i)
      end do
    end if

    print *, "############################# End Print Geometry ########################################"

    ! Write geometry to file

!    open(newunit=file_unit, file = 'raw_geometry_' // str(par_env%proc_id) // '.dat', status = 'new')  

!    write(file_unit, *) "h ", mesh%geo%h
!    write(file_unit, *) "scalefator ", mesh%geo%scalefactor
!    write(file_unit, *) "volumes ", mesh%geo%volumes(:)
!    write(file_unit, *) "face_areas ", mesh%geo%face_areas(:,:)
!    write(file_unit, *) "x_p ", mesh%geo%x_p(:,:)
!    write(file_unit, *) "x_f ", mesh%geo%x_f(:,:,:)
!    write(file_unit, *) "face_normals ", mesh%geo%face_normals(:,:,:)
!    write(file_unit, *) "vert_coords ", mesh%geo%vert_coords(:,:,:)
!    
!    close(file_unit) 

  end subroutine print_geo


  ! Print mesh topology object
  subroutine print_topo(par_env, mesh)
    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(ccs_mesh), intent(in) :: mesh       !< the mesh
    integer(ccs_int) :: i          ! loop counters
    integer(ccs_int) :: nb_elem = 10
!    integer(ccs_int) :: file_unit

    print *, "############################# Print Topology ########################################"

    print *, par_env%proc_id, "global_num_cells    : ", mesh%topo%global_num_cells
    print *, par_env%proc_id, "local_num_cells     : ", mesh%topo%local_num_cells
    print *, par_env%proc_id, "halo_num_cells      : ", mesh%topo%halo_num_cells
    print *, par_env%proc_id, "global_num_vertices : ", mesh%topo%global_num_vertices
    print *, par_env%proc_id, "vert_per_cell       : ", mesh%topo%vert_per_cell
    print *, par_env%proc_id, "total_num_cells     : ", mesh%topo%total_num_cells
    print *, par_env%proc_id, "global_num_faces    : ", mesh%topo%global_num_faces
    print *, par_env%proc_id, "num_faces           : ", mesh%topo%num_faces
    print *, par_env%proc_id, "max_faces           : ", mesh%topo%max_faces
    print *, ""

    if (allocated(mesh%topo%global_indices))    print *, par_env%proc_id, "global_indices     : ", mesh%topo%global_indices(1:nb_elem)
    if (allocated(mesh%topo%num_nb))            print *, par_env%proc_id, "num_nb             : ", mesh%topo%num_nb(1:nb_elem)
    if (allocated(mesh%topo%global_boundaries)) print *, par_env%proc_id, "global_boundaries  : ", mesh%topo%global_boundaries(1:nb_elem)
    if (allocated(mesh%topo%face_cell1))        print *, par_env%proc_id, "face_cell1         : ", mesh%topo%face_cell1(1:nb_elem)
    if (allocated(mesh%topo%face_cell2))        print *, par_env%proc_id, "face_cell2         : ", mesh%topo%face_cell2(1:nb_elem)
    if (allocated(mesh%topo%bnd_rid))           print *, par_env%proc_id, "bnd_rid            : ", mesh%topo%bnd_rid(1:nb_elem)
    if (allocated(mesh%topo%xadj))              print *, par_env%proc_id, "xadj               : ", mesh%topo%xadj(1:nb_elem)
    if (allocated(mesh%topo%adjncy))            print *, par_env%proc_id, "adjncy             : ", mesh%topo%adjncy(1:nb_elem)
    if (allocated(mesh%topo%vtxdist))           print *, par_env%proc_id, "vtxdist            : ", mesh%topo%vtxdist(1:nb_elem)
    if (allocated(mesh%topo%vwgt))              print *, par_env%proc_id, "vwgt               : ", mesh%topo%vwgt(1:nb_elem)
    if (allocated(mesh%topo%adjwgt))            print *, par_env%proc_id, "adjwgt             : ", mesh%topo%adjwgt(1:nb_elem)
    if (allocated(mesh%topo%local_partition))   print *, par_env%proc_id, "local_partition    : ", mesh%topo%local_partition(1:nb_elem)
    if (allocated(mesh%topo%global_partition))  print *, par_env%proc_id, "global_partition   : ", mesh%topo%global_partition (1:nb_elem)

    print *, ""
    if (allocated(mesh%topo%global_face_indices))   then
      do i=1, nb_elem
        print *, par_env%proc_id, "global_face_indices(1:"// str(nb_elem/2) // ", " // str(i) //")", mesh%topo%global_face_indices(1:nb_elem/2,i)
      end do
    end if

    print *, ""
    if (allocated(mesh%topo%global_vertex_indices))   then
      do i=1, nb_elem
        print *, par_env%proc_id, "global_vertex_indices(1:"// str(nb_elem/2) // ", " // str(i) //")", mesh%topo%global_vertex_indices(1:nb_elem/2,i)
      end do
    end if

    print *, ""
    if (allocated(mesh%topo%face_indices))   then
      do i=1, nb_elem
        print *, par_env%proc_id, "face_indices(1:"// str(nb_elem/2) // ", " // str(i) //")", mesh%topo%face_indices(1:nb_elem/2,i)
      end do
    end if

    print *, ""
    if (allocated(mesh%topo%nb_indices))   then
      do i=1, nb_elem
        print *, par_env%proc_id, "nb_indices(1:"// str(nb_elem/2) // ", " // str(i) //")", mesh%topo%nb_indices(1:nb_elem/2,i)
      end do
    end if

    print *, "############################# End Print Topology ########################################"

!
!    open(newunit=file_unit, file = 'raw_topology_' // str(par_env%proc_id) // '.dat', status = 'new')  
!
!    write(file_unit, *) "global_num_cells         ", mesh%topo%global_num_cells
!    write(file_unit, *) "local_num_cells          ", mesh%topo%local_num_cells
!    write(file_unit, *) "halo_num_cells           ", mesh%topo%halo_num_cells
!    write(file_unit, *) "global_num_vertices      ", mesh%topo%global_num_vertices
!    write(file_unit, *) "vert_per_cell            ", mesh%topo%vert_per_cell
!    write(file_unit, *) "total_num_cells          ", mesh%topo%total_num_cells
!    write(file_unit, *) "global_num_faces         ", mesh%topo%global_num_faces
!    write(file_unit, *) "num_faces                ", mesh%topo%num_faces
!    write(file_unit, *) "max_faces                ", mesh%topo%max_faces
!    write(file_unit, *) "global_indices           ", mesh%topo%global_indices
!    write(file_unit, *) "global_face_indices      ", mesh%topo%global_face_indices
!    write(file_unit, *) "global_vertex_indices    ", mesh%topo%global_vertex_indices
!    write(file_unit, *) "face_indices             ", mesh%topo%face_indices
!    write(file_unit, *) "nb_indices               ", mesh%topo%nb_indices
!    write(file_unit, *) "num_nb                   ", mesh%topo%num_nb
!    write(file_unit, *) "global_boundaries        ", mesh%topo%global_boundaries
!    write(file_unit, *) "face_cell1               ", mesh%topo%face_cell1
!    write(file_unit, *) "face_cell2               ", mesh%topo%face_cell2
!    write(file_unit, *) "bnd_rid                  ", mesh%topo%bnd_rid
!    write(file_unit, *) "xadj                     ", mesh%topo%xadj
!    write(file_unit, *) "adjncy                   ", mesh%topo%adjncy
!    write(file_unit, *) "vtxdist                  ", mesh%topo%vtxdist
!    write(file_unit, *) "vwgt                     ", mesh%topo%vwgt
!    write(file_unit, *) "adjwgt                   ", mesh%topo%adjwgt
!    write(file_unit, *) "local_partition          ", mesh%topo%local_partition
!    write(file_unit, *) "global_partition         ", mesh%topo%global_partition     
!
!    close(file_unit)
!
  end subroutine print_topo


end module mesh_utils
