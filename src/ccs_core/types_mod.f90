!v Module file types.mod
!
!  Provides concrete types and bases of extensible types.

module types

  use kinds, only: ccs_int, ccs_real, ccs_long
  use parallel_types, only: parallel_environment

  implicit none

  private

  !> Stub type for vectors to be extended in sub-modules.
  type, public :: ccs_vector
    character(len=:), allocatable :: name  !< Name of the vector object
  end type ccs_vector

  !> Stub type for matrices to be extended in sub-modules.
  type, public :: ccs_matrix
    character(len=:), allocatable :: name  !< Name of the matrix object
  end type ccs_matrix

  !> Container type for data required to initialise a vector.
  type, public :: vector_spec
    class(parallel_environment), pointer :: par_env !< The parallel environment
    type(ccs_mesh), pointer :: mesh                 !< The mesh object to build the vector on
    integer(ccs_int) :: storage_location            !< The storage location of the vector values (cell or face)
  end type vector_spec

  !> Container type for setting values in a vector.
  type, public :: vector_values
    integer(ccs_int), dimension(:), allocatable :: global_indices !< Array of (global) indices to set values
    !< on, must be same size as values array.
    real(ccs_real), dimension(:), allocatable :: values           !< Array of values, must be same size as
    !< index array.
    integer(ccs_int) :: setter_mode                               !< Which mode to use when setting values?
    integer(ccs_int) :: current_entry                             !< Which entry are we currently working on?
  end type vector_values

  !> Container type for data required to initialise a matrix.
  type, public :: matrix_spec
    type(ccs_mesh), pointer :: mesh                 !< The mesh
    class(parallel_environment), pointer :: par_env !< The parallel environment
    integer(ccs_int) :: nnz                         !< Non-zeros per row
  end type matrix_spec

  !> Container type for setting values in a matrix.
  type, public :: matrix_values
    integer(ccs_int), dimension(:), allocatable :: global_row_indices !< Array of (global) row indices to set values on.
    integer(ccs_int), dimension(:), allocatable :: global_col_indices !< Array of (global) column indices to set values on.
    real(ccs_real), dimension(:), allocatable :: values     !< Array of values, must be logically 2D and
    !< of size = size(row_indices) * size(col_indices). Uses
    !< row-major ordering.
    integer(ccs_int) :: setter_mode                         !< Which mode to use when setting values?
    integer(ccs_int) :: current_row, current_col            !< Which entry are we currently working on?
  end type matrix_values

  type, public :: matrix_values_spec
    integer(ccs_int) :: nrows = 0
    integer(ccs_int) :: ncols = 0
  end type matrix_values_spec

  !>  Container type representing a linear system.
  type, public :: equation_system
    character(len=:), allocatable :: name  !< Name of the equation system
    class(ccs_vector), pointer :: solution !< Solution vector
    class(ccs_vector), pointer :: rhs      !< Right-hand side vector
    class(ccs_matrix), pointer :: matrix   !< Matrix
    class(parallel_environment), pointer :: par_env !< The parallel environment
  end type equation_system

  !> Stub type for solvers to be extended in sub-modules.
  type, public :: linear_solver
    type(equation_system) :: linear_system !< System of equations
  end type linear_solver

  !v Graph connectivity type
  type, public :: graph_connectivity
    integer(ccs_long), dimension(:), allocatable :: xadj            !< Array that points to where in adjncy the list for each vertex
                                                                    !<   begins and ends  - name from ParMETIS
    integer(ccs_long), dimension(:), allocatable :: adjncy          !< Array storing adjacency lists for each vertex consecutively
                                                                    !<   - name from ParMETIS
    integer(ccs_long), dimension(:), allocatable :: vtxdist         !< Array that indicates vertices local to a processor. Rank p_i stores
                                                                    !<   the vertices from vtxdist[i] up to (but not including) vertex
                                                                    !<   vtxdist[i + 1] - name from ParMETIS
    integer(ccs_long), dimension(:), allocatable :: vwgt            !< Weights on vertices - name from ParMETIS
    integer(ccs_long), dimension(:), allocatable :: adjwgt          !< Weights on edges - name from ParMETIS
    integer(ccs_long), dimension(:), allocatable :: local_partition !< Local partition array
    integer(ccs_long), dimension(:), pointer :: global_partition    !< Global partition array
    integer :: global_partition_window                              !< Associated shared window
  end type graph_connectivity
  
  !v Topology type
  !
  !  Describes the topology (i.e. connectivity) of the mesh.
  !  This includes the numbering of cells, for which there are 3 important values:
  !  - local index: this is implicitly defined by 1 <= local_index <= n, with halo cells stored at local indices > local_num_cells,
  !  - natural index: this is the global index of the cells as originally defined by the mesh, a given process may have
  !                   discontiguous range(s) of natural indices.
  !  - global index: the index (i.e. row) of a cell in the linear system. Each process has a contiguous range of global indices,
  !                  i.e. each global index is given by the local index + a constant pre-process offset.
  type, public :: topology
    type(graph_connectivity) :: graph_conn                                  !< Object describing the connectivity
    integer(ccs_int) :: global_num_cells                                    !< Global number of cells
    integer(ccs_int) :: local_num_cells                                     !< Local number of cells
    integer(ccs_int) :: halo_num_cells                                      !< Local number of halo cells
    integer(ccs_int) :: global_num_vertices                                 !< Global number of vertices
    integer(ccs_int) :: vert_per_cell                                       !< Number of vertices per cell
    integer(ccs_int) :: total_num_cells                                     !< Number of local + halo cells
    integer(ccs_int) :: global_num_faces                                    !< Global number of faces
    integer(ccs_int) :: num_faces                                           !< Local number of faces
    integer(ccs_int) :: max_faces                                           !< Maximum number of faces per cell
    integer(ccs_int), dimension(:), allocatable :: natural_indices          !< The global index of cells in the original ordering (local + halo)
                                                                            !<   natural_icell = natural_indices(local_icell)
    integer(ccs_int), dimension(:), allocatable :: global_indices           !< The global index of cells (local + halo)
                                                                            !<   global_icell = global_indices(local_icell)
    integer(ccs_int), dimension(:, :), pointer :: global_face_indices       !< Global list of faces indices
                                                                            !<   global_iface = global_face_indices(cell_iface, global_icell)
                                                                            !<   (no special treatment for halo or boundary faces)
    integer :: global_face_indices_window                                   !< Associated shared window
    integer(ccs_int), dimension(:, :), allocatable :: loc_global_vertex_indices     !< local version of the global list of vertex indices
                                                                            !<   global_ivert = loc_global_vertex_indices(ivert, local_icell)
    integer(ccs_int), dimension(:, :), pointer :: global_vertex_indices     !< Global list of vertex indices
                                                                            !<   global_ivert = global_vertex_indices(ivert, global_icell)
    integer :: global_vertex_indices_window                                 !< Associated shared window
    integer(ccs_int), dimension(:, :), allocatable :: face_indices          !< Cell face index in local face vector (face, cell)
                                                                            !<   iface = global_face_indices(cell_iface, icell)
                                                                            !<   (no special treatment for halo or boundary faces)
    integer(ccs_int), dimension(:, :), allocatable :: nb_indices            !< Cell face index in local face vector (face, cell)
                                                                            !<   nb_icell = nb_indices(cell_iface, icell) -> returns <0 on boundaries
    integer(ccs_int), dimension(:), allocatable :: num_nb                   !< The local number of neighbours per cell
                                                                            !<   num_nb = num_nb(icell), equiv to number of faces, boundary 'neighbours' are counted
    integer(ccs_int), dimension(:), pointer :: face_cell1                   !< Array of 1st face cells
                                                                            !<   global_icell1 = face_cell1(global_iface).
    integer :: face_cell1_window                                            !< Associated shared window
    integer(ccs_int), dimension(:), pointer :: face_cell2                   !< Array of 2nd face cells
                                                                            !<   global_icell2 = face_cell2(global_iface) -> returns 0 on boundaries
    integer :: face_cell2_window                                            !< Associated shared window
    integer(ccs_int), dimension(:), pointer :: bnd_rid                      !< global face boundary index.
                                                                            !< 0 on internal faces
                                                                            !< -X on a bondary face according to the boundary index
    integer :: bnd_rid_window                                               !< Associated shared window
    integer(ccs_int) :: shared_array_local_offset                           !< Offset within shared arrays for quantities that are locally indexed (i.e. each rank is responsible for local_num_cells of these)
    integer(ccs_int) :: shared_array_total_offset                           !< Offset within shared arrays for quantities that are totally indexed (i.e. each rank is responsible for total_num_cells of these)
  end type topology

  !> Geometry type
  type, public :: geometry
    real(ccs_real) :: h                                                 !< The (constant) grid spacing XXX: remove!
    real(ccs_real) :: scalefactor                                       !< Scalefactor
    real(ccs_real), dimension(:, :), pointer :: face_areas              !< Face areas (face, cell)
    integer :: face_areas_window                                        !< Associated shared window
    real(ccs_real), dimension(:), pointer :: volumes                    !< Cell volumes
    integer :: volumes_window                                           !< Associated shared window
    real(ccs_real), dimension(:, :), pointer :: x_p                     !< Cell centres (dimension, cell)
    integer :: x_p_window                                               !< Associated shared window
    real(ccs_real), dimension(:, :, :), pointer :: x_f                  !< Face centres (dimension, face, cell)
    integer :: x_f_window                                               !< Associated shared window
    real(ccs_real), dimension(:, :, :), pointer :: face_normals         !< Face normals (dimension, face, cell)
    integer :: face_normals_window                                      !< Associated shared window
    real(ccs_real), dimension(:, :, :), pointer :: vert_coords          !< Vertex coordinates (dimension, vertex, cell)
    integer :: vert_coords_window                                       !< Associated shared window
    real(ccs_real), dimension(:), allocatable :: face_interpol          !< Face interpolation factor, factor = face_interpol(iface)
  end type geometry

  !> Mesh type
  type, public :: ccs_mesh
    type(topology) :: topo
    type(geometry) :: geo
    logical :: is_generated                                    !< Indicates whether mesh was generated (true) or read (false)
    character(len=128), dimension(:), allocatable :: bnd_names !< Array of boundary names, index corresponding to the boundary ID
  end type ccs_mesh

  !> BC data type
  type, public :: bc_config
    integer(ccs_int), dimension(:), allocatable :: ids
    integer(ccs_int), dimension(:), allocatable :: bc_types
    real(ccs_real), dimension(:), allocatable :: values
    type(bc_profile), dimension(:), allocatable :: profiles
  end type bc_config

  !> Boundary condition profile
  type, public :: bc_profile
    real(ccs_real), dimension(:), allocatable :: centre      ! reference location for the coordinates
    real(ccs_real), dimension(:), allocatable :: coordinates ! coordinates at which a value is defined, this is the distance to the centre reference location
    real(ccs_real), dimension(:), allocatable :: values      ! boundary condition values associated to each coordinate
  end type bc_profile

  !v Wrapper class for ccs_vector
  !
  !  A wrapper is required for ccs_vector in order to allow the creation of
  !  an allocatable array of ccs_vectors. ccs_vectors need to be allocatable themselves,
  !  but it is not possible to have an allocatable array with allocatable elements. Having a wrapper
  !  gets around this issue by essentially enabling the creation of a dynamic array of pointers
  !  which point to ccs_vector elements. This is used in, e.g. the fields datatype to store
  !  (an unknown number of) past values for use by timestepping schemes.
  type :: ccs_vector_ptr
    class(ccs_vector), allocatable :: vec
  end type ccs_vector_ptr

  !> Scalar field type
  type, public :: field
    class(ccs_vector), allocatable :: values                      !< Vector representing the field
    class(ccs_vector), allocatable :: residuals                   !< Vector representing the field's residuals
    type(ccs_vector_ptr), dimension(:), allocatable :: old_values !< Vector representing the old fields
    class(ccs_vector), allocatable :: x_gradients                 !< Vector representing the x gradient
    class(ccs_vector), allocatable :: y_gradients                 !< Vector representing the y gradient
    class(ccs_vector), allocatable :: z_gradients                 !< Vector representing the z gradient
    real(ccs_real), dimension(:), pointer :: values_ro            !< Read only pointer to array containing values
    real(ccs_real), dimension(:), pointer :: x_gradients_ro       !< Read only pointer to array containing x_gradients
    real(ccs_real), dimension(:), pointer :: y_gradients_ro       !< Read only pointer to array containing y_gradients
    real(ccs_real), dimension(:), pointer :: z_gradients_ro       !< Read only pointer to array containing z_gradients
    type(bc_config) :: bcs                                        !< The bcs data structure for the cell
    real(ccs_real) :: Schmidt = 1.0                               !< Schmidt Number
    logical :: enable_cell_corrections                            !< Whether or not deffered corrections should be used (non-orthogonality, excentricity etc.)
    character(len=20) :: name
    logical :: output = .false.                                   !< Should field be written in output?
    logical :: solve = .true.                                     !< Whether to solve a linear system for this variable or not
  end type field

  type, public, extends(field) :: upwind_field
  end type
  type, public, extends(field) :: central_field
  end type
  type, public, extends(field) :: face_field
  end type
  type, public, extends(field) :: gamma_field
  end type
  type, public, extends(field) :: linear_upwind_field
  end type

  !> Field specification type, used for defining new fields.
  type, public :: field_spec
    character(len=:), allocatable :: ccs_config_file !< Config file containing field information
    type(vector_spec) :: vec_properties              !< Descriptor for the underlying vector
    integer :: field_type                            !< Flag to identify which type of field to create
    character(len=:), allocatable :: field_name      !< The name of the field
    integer(ccs_int) :: n_boundaries                 !< The number of boundaries involved...
    logical :: store_residuals = .false.             !< Whether or not residuals should be stored for this field
    logical :: enable_cell_corrections = .true.      !< Whether or not deffered corrections should be used (non-orthogonality, excentricity etc.)
  end type field_spec

  !> Type for storing pointer to a field
  type, public :: field_ptr
    class(field), pointer :: ptr => null()   !< Pointer to the field data
    character(len=:), allocatable :: name    !< Name of the field
  end type field_ptr

  !v Cell locator
  !
  ! Lightweight type to provide easy cell location based on a cell's cell connectivity.
  type, public :: cell_locator
    integer(ccs_int) :: index_p     !< Cell index
  end type cell_locator

  !v Face locator
  !
  !  Lightweight type to provide easy face location based on a cell's face connectivity.
  type, public :: face_locator
    integer(ccs_int) :: index_p       !< Cell index
    integer(ccs_int) :: cell_face_ctr !< Cell-face ctr i.e. I want to access face "3" of the cell.
  end type face_locator

  !v Neighbour locator
  !
  !  Lightweight type to provide easy cell-neighbour connection.
  type, public :: neighbour_locator
    integer(ccs_int) :: index_p       !< the cell index relative to which this is a neighbour
    integer(ccs_int) :: nb_counter    !< the cell-relative counter identifying this neighbour
  end type neighbour_locator

  !v Vertex locator
  !
  !  Lightweight type to provide easy vertex location based on a cell's vertex connectivity.
  type, public :: vert_locator
    integer(ccs_int) :: index_p       !< Cell index
    integer(ccs_int) :: cell_vert_ctr !< Cell-vertex ctr i.e. I want to access vertex "3" of the cell.
  end type vert_locator

  !v Fluid type
  !
  ! Type for accumulating all the fluid data
  type, public :: fluid
    type(field_ptr), dimension(:), allocatable :: fields
  end type fluid

  !>  IO environment type
  type, public :: io_environment
  end type io_environment

  !> Process that will perform file IO
  type, public :: io_process
  end type io_process

end module types
