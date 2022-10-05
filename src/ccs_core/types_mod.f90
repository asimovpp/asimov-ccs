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
  end type ccs_vector

  !> Stub type for matrices to be extended in sub-modules.
  type, public :: ccs_matrix
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

  !> Topology type
  type, public :: topology
    integer(ccs_int) :: global_num_cells                              !< Global number of cells
    integer(ccs_int) :: local_num_cells                               !< Local number of cells
    integer(ccs_int) :: halo_num_cells                                !< Number of halo cells
    integer(ccs_int) :: global_num_vertices                           !< Global number of vertices
    integer(ccs_int) :: vert_per_cell                                 !< Number of vertices per cell
    integer(ccs_int) :: total_num_cells                               !< Number of local + halo cells        
    integer(ccs_int) :: global_num_faces                              !< Global number of faces
    integer(ccs_int) :: num_faces                                     !< Local number of faces
    integer(ccs_int) :: max_faces                                     !< Maximum number of faces per cell
    integer(ccs_int), dimension(:), allocatable :: global_indices         !< The global index of cells (local + halo)
    integer(ccs_int), dimension(:, :), allocatable :: global_face_indices !< Global list of faces indices
    integer(ccs_int), dimension(:, :), allocatable :: global_vertex_indices !< Global list of vertex indices
    integer(ccs_int), dimension(:, :), allocatable :: face_indices        !< Cell face index in local face vector (face, cell)
    integer(ccs_int), dimension(:, :), allocatable :: nb_indices      !< Cell face index in local face vector (face, cell)
    integer(ccs_int), dimension(:), allocatable :: num_nb             !< The local number of neighbours per cell
    integer(ccs_int), dimension(:), allocatable :: global_boundaries  !< Array of boundary faces
    integer(ccs_int), dimension(:), allocatable :: face_cell1         !< Array of 1st face cells
    integer(ccs_int), dimension(:), allocatable :: face_cell2         !< Array of 2nd face cells
    integer(ccs_long), dimension(:), allocatable :: xadj              !< Array that points to where in adjncy the list for each vertex
    !< begins and ends  - name from ParMETIS
    integer(ccs_long), dimension(:), allocatable :: adjncy            !< Array storing adjacency lists for each vertex consecutively
    !< - name from ParMETIS
    integer(ccs_long), dimension(:), allocatable :: vtxdist           !< Array that indicates vertices local to a processor. Rank p_i stores
    !< the vertices from vtxdist[i] up to (but not including) vertex
    !< vtxdist[i + 1] - name from ParMETIS
    integer(ccs_long), dimension(:), allocatable :: vwgt              !< Weights on vertices - name from ParMETIS
    integer(ccs_long), dimension(:), allocatable :: adjwgt            !< Weights on edges - name from ParMETIS
    integer(ccs_long), dimension(:), allocatable :: local_partition   !< Local partition array
    integer(ccs_long), dimension(:), allocatable :: global_partition  !< Global partition array
  end type topology

  !> Geometry type
  type, public :: geometry
    real(ccs_real) :: h                                                 !< The (constant) grid spacing XXX: remove!
    real(ccs_real) :: scalefactor                                       !< Scalefactor 
    real(ccs_real), dimension(:, :), allocatable :: face_areas          !< Face areas (face, cell)
    real(ccs_real), dimension(:), allocatable :: volumes                !< Cell volumes
    real(ccs_real), dimension(:, :), allocatable :: x_p                 !< Cell centres (dimension, cell)
    real(ccs_real), dimension(:, :, :), allocatable :: x_f              !< Face centres (dimension, face, cell)
    real(ccs_real), dimension(:, :, :), allocatable :: face_normals     !< Face normals (dimension, face, cell)
    real(ccs_real), dimension(:, :, :), allocatable :: vert_coords      !< Vertex coordinates (dimension, vertex, cell)
  end type geometry

  !> Mesh type
  type, public :: ccs_mesh
    type(topology) :: topo
    type(geometry) :: geo
  end type ccs_mesh

  !> BC data type
  type, public :: bc_config
    integer(ccs_int), dimension(:), allocatable :: ids
    integer(ccs_int), dimension(:), allocatable :: bc_types
    real(ccs_real), dimension(:), allocatable :: values
  end type bc_config

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
    type(ccs_vector_ptr), dimension(:), allocatable :: old_values !< Vector representing the old fields
    class(ccs_vector), allocatable :: x_gradients                 !< Vector representing the x gradient
    class(ccs_vector), allocatable :: y_gradients                 !< Vector representing the y gradient
    class(ccs_vector), allocatable :: z_gradients                 !< Vector representing the z gradient
    type(bc_config) :: bcs                                        !< The bcs data structure for the cell
  end type field

  type, public, extends(field) :: upwind_field
  end type
  type, public, extends(field) :: central_field
  end type
  type, public, extends(field) :: face_field
  end type

  !v Cell locator
  !
  ! Lightweight type to provide easy cell location based on a cell's cell connectivity.
  type, public :: cell_locator
    type(ccs_mesh), pointer :: mesh !< Pointer to the mesh -- we DON'T want to copy this!
    integer(ccs_int) :: index_p     !< Cell index
  end type cell_locator

  !v Face locator
  !
  !  Lightweight type to provide easy face location based on a cell's face connectivity.
  type, public :: face_locator
    type(ccs_mesh), pointer :: mesh   !< Pointer to the mesh -- we DON'T want to copy this!
    integer(ccs_int) :: index_p       !< Cell index
    integer(ccs_int) :: cell_face_ctr !< Cell-face ctr i.e. I want to access face "3" of the cell.
  end type face_locator

  !v Neighbour locator
  !
  !  Lightweight type to provide easy cell-neighbour connection.
  type, public :: neighbour_locator
    type(ccs_mesh), pointer :: mesh
    integer(ccs_int) :: index_p
    integer(ccs_int) :: nb_counter
  end type neighbour_locator

  !>  IO environment type
  type, public :: io_environment
  end type io_environment

  !> Process that will perform file IO
  type, public :: io_process
  end type io_process

end module types
