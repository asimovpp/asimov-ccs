!v Module file types.mod
!
!  Provides concrete types and bases of extensible types.

module types

  use kinds, only: ccs_int, ccs_real
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

  !> BC data type
  type, public :: bc_config
    integer(ccs_int), dimension(:), allocatable :: ids
    integer(ccs_int), dimension(:), allocatable :: bc_types
    real(ccs_real), dimension(:), allocatable :: values
  end type bc_config

  !> Mesh type
  type, public :: ccs_mesh
    integer(ccs_int) :: nglobal      !< Global mesh size
    integer(ccs_int) :: nlocal       !< Local mesh size
    integer(ccs_int) :: nhalo        !< How many cells in my halo?
    integer(ccs_int) :: ntotal       !< How many cells do I interact with (nlocal + nhalo)?
    integer(ccs_int) :: nfaces_local !< Number of faces in local mesh
    integer(ccs_int), dimension(:), allocatable :: global_indices       !< The global index of cells (local + halo)
    integer(ccs_int), dimension(:), allocatable :: nnb                  !< The per-cell neighbour count
    integer(ccs_int), dimension(:, :), allocatable :: neighbour_indices !< Cell neighbours (neighbour/face, cell)
    integer(ccs_int), dimension(:, :), allocatable :: face_indices      !< Cell face index in local face vector (face, cell)
    real(ccs_real) :: h                                                 !< The (constant) grid spacing XXX: remove!
    real(ccs_real), dimension(:, :), allocatable :: face_areas          !< Face areas
    real(ccs_real), dimension(:), allocatable :: volumes                !< Cell volumes
    real(ccs_real), dimension(:, :), allocatable :: x_p                 !< Cell centres (dimension, cell)
    real(ccs_real), dimension(:, :, :), allocatable :: x_f              !< Face centres (dimension, face, cell)
    real(ccs_real), dimension(:, :, :), allocatable :: face_normals     !< Face normals (dimension, face, cell)
  end type ccs_mesh

  !> Scalar field type
  type, public :: field
    class(ccs_vector), allocatable :: values      !< Vector representing the field
    class(ccs_vector), allocatable :: x_gradients !< Vector representing the x gradient
    class(ccs_vector), allocatable :: y_gradients !< Vector representing the y gradient
    class(ccs_vector), allocatable :: z_gradients !< Vector representing the z gradient
    type(bc_config) :: bcs                        !< The bcs data structure for the cell
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
