!> @brief Module file types.mod
!
!> @details Provides concrete types and bases of extensible types.

module types

  use kinds, only : ccs_int, ccs_real
  use parallel_types, only: parallel_environment
  
  implicit none

  private

  !> @brief Stub type for vectors to be extended in sub-modules.
  type, public :: ccs_vector
  end type ccs_vector

  !> @brief Stub type for matrices to be extended in sub-modules.
  type, public :: ccs_matrix
  end type ccs_matrix

  !> @brief Container type for data required to initialise a vector.
  type, public :: vector_spec
    class(parallel_environment), pointer :: par_env !< The parallel environment
    type(ccs_mesh), pointer :: mesh !< The mesh object to build the vector on
    integer(ccs_int) :: storage_location !< The storage location of the vector values (cell or face)
  end type vector_spec


  !> @brief Container type for setting values in a vector.
  type, public :: vector_values
    integer(ccs_int), dimension(:), allocatable :: indices !< Array of (global) indices to set values
                                                        !! on, must be same size as values array.
    real(ccs_real), dimension(:), allocatable :: values   !< Array of values, must be same size as
                                                        !! index array.
    integer(ccs_int) :: setter_mode                           !< Which mode to use when setting values?
  end type vector_values

  !> @brief Container type for data required to initialise a matrix.
  type, public :: matrix_spec
    type(ccs_mesh), pointer :: mesh                     !< The mesh
    class(parallel_environment), pointer :: par_env !< The parallel environment
    integer(ccs_int) :: nnz                        !< Non-zeros per row
  end type matrix_spec

  !> @brief Container type for setting values in a matrix.
  type, public :: matrix_values
    integer(ccs_int), dimension(:), allocatable :: row_indices !< Array of (global) row indices to set values on.
    integer(ccs_int), dimension(:), allocatable :: col_indices !< Array of (global) column indices to set values on.
    real(ccs_real), dimension(:), allocatable :: values     !< Array of values, must be logically 2D and 
                                                          !! of size = size(row_indices) * size(col_indices). Uses 
                                                          !! row-major ordering.
    integer(ccs_int) :: setter_mode !< Which mode to use when setting values?
  end type matrix_values

  !> @brief Container type representing a linear system.
  type, public :: equation_system
    class(ccs_vector), pointer :: solution !< Solution vector
    class(ccs_vector), pointer :: rhs !< Right-hand side vector
    class(ccs_matrix), pointer :: matrix   !< Matrix
    class(parallel_environment), pointer :: par_env !< The parallel environment
  end type equation_system
  
  !> @brief Stub type for solvers to be extended in sub-modules.
  type, public :: linear_solver
    type(equation_system) :: linear_system !< System of equations
  end type linear_solver

  !> @brief Mesh type
  type, public :: ccs_mesh
    integer(ccs_int) :: nglobal !< Global mesh size
    integer(ccs_int) :: nlocal  !< Local mesh size
    integer(ccs_int) :: nhalo   !< How many cells in my halo?
    integer(ccs_int) :: ntotal  !< How many cells do I interact with (nlocal + nhalo)?
    integer(ccs_int) :: nfaces_local !< Number of faces in local mesh
    integer(ccs_int), dimension(:), allocatable :: idx_global ! The global index of cells (local + halo)
    integer(ccs_int), dimension(:), allocatable :: nnb        ! The per-cell neighbour count
    integer(ccs_int), dimension(:, :), allocatable :: index_nb !< Cell neighbours (neighbour/face, cell)
    integer(ccs_int), dimension(:, :), allocatable :: faceidx  !< Cell face index in local face vector (face, cell)
    real(ccs_real) :: h                                        !< The (constant) grid spacing XXX: remove!
    real(ccs_real), dimension(:, :), allocatable :: Af         !< Face areas
    real(ccs_real), dimension(:), allocatable :: vol           !< Cell volumes
    real(ccs_real), dimension(:, :), allocatable :: xc         !< Cell centres (dimension, cell)
    real(ccs_real), dimension(:, :, :), allocatable :: xf      !< Face centres (dimension, face, cell)
    real(ccs_real), dimension(:, :, :), allocatable :: nf      !< Face normals (dimension, face, cell)
  end type ccs_mesh

  !> @brief Scalar field type
  type, public :: field
    class(ccs_vector), allocatable :: values   !< Vector representing the field
    class(ccs_vector), allocatable :: x_gradients !< Vector representing the x gradient
    class(ccs_vector), allocatable :: y_gradients !< Vector representing the y gradient
    class(ccs_vector), allocatable :: z_gradients !< Vector representing the z gradient
  end type field

  type, public, extends(field) :: upwind_field
  end type
  type, public, extends(field) :: central_field
  end type
  type, public, extends(field) :: face_field
  end type

  !> @brief Cell locator
  !
  !> @description Lightweight type to provide easy cell location based on a cell's cell
  !!              connectivity.
  type, public :: cell_locator
    type(ccs_mesh), pointer :: mesh        !< Pointer to the mesh -- we DON'T want to copy this!
    integer(ccs_int) :: index_p      !< Cell index
  end type cell_locator

  !> @brief Face locator
  !
  !> @description Lightweight type to provide easy face location based on a cell's face
  !!              connectivity.
  type, public :: face_locator
    type(ccs_mesh), pointer :: mesh        !< Pointer to the mesh -- we DON'T want to copy this!
    integer(ccs_int) :: index_p      !< Cell index
    integer(ccs_int) :: cell_face_ctr !< Cell-face ctr i.e. I want to access face "3" of the cell.
  end type face_locator

  !> @brief Neighbour locator
  !
  !> @description Lightweight type to provide easy cell-neighbour connection.
  type, public :: neighbour_locator
    type(ccs_mesh), pointer :: mesh
    integer(ccs_int) :: index_p
    integer(ccs_int) :: nb_counter
  end type neighbour_locator

  type, public :: bc_config
    integer(ccs_int), dimension(4) :: region
    integer(ccs_int), dimension(4) :: bc_type
    real(ccs_real), dimension(4, 2) :: endpoints ! Used in scalar_advection case and tests, 
                                                  ! possibly remove/improve for general
  end type bc_config
  
  !> @brief IO environment type
  type, public :: io_environment
  end type io_environment

  !> @brief Process that will perform file IO
  type, public :: io_process
  end type io_process

end module types
