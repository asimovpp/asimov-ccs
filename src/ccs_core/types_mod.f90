!> @brief Module file types.mod
!
!> @details Provides concrete types and bases of extensible types.

module types

  use kinds, only : accs_int, accs_real
  use parallel_types, only: parallel_environment
  
  implicit none

  private

  public :: set_global_matrix_size
  public :: set_local_matrix_size
  public :: set_face_location
  public :: set_cell_location
  public :: set_neighbour_location

  type, public :: viewer
  end type viewer

  !> @brief Stub type for vectors to be extended in sub-modules.
  type, public :: vector
  end type vector

  !> @brief Stub type for matrices to be extended in sub-modules.
  type, public :: matrix
  end type matrix

  !> @brief Container type for data required to initialise a vector.
  type, public :: vector_init_data
    !> The vector size can be specified either globally or per-process
    integer(accs_int) :: nglob !> The global vector size (set -1 to ignore)
    integer(accs_int) :: nloc  !> The local vector size (set -1 to ignore)
    class(parallel_environment), pointer :: par_env !> The parallel environment
  end type vector_init_data


  !> @brief Container type for setting values in a vector.
  type, public :: vector_values
    integer(accs_int), dimension(:), allocatable :: idx !> Array of (global) indices to set values
                                                        !! on, must be same size as values array.
    real(accs_real), dimension(:), allocatable :: val   !> Array of values, must be same size as
                                                        !! index array.
    integer(accs_int) :: mode                           !> Which mode to use when setting values?
  end type vector_values

  !> @brief Container type for data required to initialise a matrix.
  type, public :: matrix_init_data
    !> The matrix size can be specified either globally or per-process
    integer(accs_int) :: rglob !> The global matrix rows size (set -1 to ignore)
    integer(accs_int) :: cglob !> The global matrix columns size (set -1 to ignore)
    integer(accs_int) :: rloc  !> The local matrix rows size (set -1 to ignore)
    integer(accs_int) :: cloc  !> The local matrix columns size (set -1 to ignore)
    integer(accs_int) :: nnz   !> The number of non-zeros in a row - MUST include the 
                              !! diagonal regardles off value.
    class(parallel_environment), pointer :: par_env !> The parallel environment
  end type matrix_init_data

  !> @brief Container type for setting values in a matrix.
  type, public :: matrix_values
    integer(accs_int), dimension(:), allocatable :: rglob !> Array of (global) row indices to set values on.
    integer(accs_int), dimension(:), allocatable :: cglob !> Array of (global) column indices to set values on.
    real(accs_real), dimension(:), allocatable :: val     !> Array of values, must be logically 2D and 
                                                          !! of size = size(rglob) * size(cglob). Uses 
                                                          !! row-major ordering.
    integer(accs_int) :: mode !> Which mode to use when setting values?
  end type matrix_values

  !> @brief Container type representing a linear system.
  type, public :: linear_system
    class(vector), pointer :: sol !> Solution vector
    class(vector), pointer :: rhs !> Right-hand side vector
    class(matrix), pointer :: M   !> Matrix
    class(parallel_environment), pointer :: par_env !> The parallel environment
  end type linear_system
  
  !> @brief Stub type for solvers to be extended in sub-modules.
  type, public :: linear_solver
    type(linear_system) :: eqsys !> System of equations
  end type linear_solver

  !> @brief Mesh type
  type, public :: mesh
    integer(accs_int) :: n !> Global mesh size
    integer(accs_int) :: nlocal !> Local mesh size
    integer(accs_int), dimension(:), allocatable :: idx_global 
    integer(accs_int), dimension(:), allocatable :: nnb 
    integer(accs_int), dimension(:, :), allocatable :: nbidx !> Cell neighbours (neighbour/face, cell)
    real(accs_real) :: h
    real(accs_real), dimension(:, :), allocatable :: Af      !> Face areas
    real(accs_real), dimension(:), allocatable :: vol        !> Cell volumes
    real(accs_real), dimension(:, :), allocatable :: xc      !> Cell centres (dimension, cell)
    real(accs_real), dimension(:, :, :), allocatable :: xf   !> Face centres (dimension, face, cell)
    real(accs_real), dimension(:, :, :), allocatable :: nf   !> Face normals (dimension, face, cell)
  end type mesh

  !> @brief Scalar field type
  type, public :: field
    real(accs_real), dimension(:,:), allocatable :: val
  end type field

  type, public, extends(field) :: upwind_field
  end type
  type, public, extends(field) :: central_field
  end type

  !> @brief Cell locator
  !
  !> @description Lightweight type to provide easy cell location based on a cell's cell
  !!              connectivity.
  type, public :: cell_locator
    type(mesh), pointer :: mesh        !> Pointer to the mesh -- we DON'T want to copy this!
    integer(accs_int) :: cell_idx      !> Cell index
  end type cell_locator

  !> @brief Face locator
  !
  !> @description Lightweight type to provide easy face location based on a cell's face
  !!              connectivity.
  type, public :: face_locator
    type(mesh), pointer :: mesh        !> Pointer to the mesh -- we DON'T want to copy this!
    integer(accs_int) :: cell_idx      !> Cell index
    integer(accs_int) :: cell_face_ctr !> Cell-face ctr i.e. I want to access face "3" of the cell.
  end type face_locator

  !> @brief Neighbour locator
  !
  !> @description Lightweight type to provide easy cell-neighbour connection.
  type, public :: neighbour_locator
    type(mesh), pointer :: mesh
    integer(accs_int) :: cell_idx
    integer(accs_int) :: cell_neighbour_ctr
  end type neighbour_locator
  
  interface
  module subroutine set_global_matrix_size(mat, rows, columns, nnz, par_env)
    type(matrix_init_data), intent(inout) :: mat
    integer(accs_int), intent(in) :: rows
    integer(accs_int), intent(in) :: columns
    integer(accs_int), intent(in) :: nnz
    class(parallel_environment), allocatable, target, intent(in) :: par_env
  end subroutine set_global_matrix_size

  module subroutine set_local_matrix_size(mat, rows, columns, nnz, par_env)
    type(matrix_init_data), intent(inout) :: mat
    integer(accs_int), intent(in) :: rows
    integer(accs_int), intent(in) :: columns
    integer(accs_int), intent(in) :: nnz
    class(parallel_environment), allocatable, target, intent(in) :: par_env
  end subroutine set_local_matrix_size

  module subroutine set_cell_location(cell_location, geometry, cell_idx)
    type(cell_locator), intent(out) :: cell_location
    type(mesh), target, intent(in) :: geometry
    integer(accs_int), intent(in) :: cell_idx
  end subroutine set_cell_location

  module subroutine set_face_location(face_location, geometry, cell_idx, cell_face_ctr)
    type(face_locator), intent(out) :: face_location
    type(mesh), target, intent(in) :: geometry
    integer(accs_int), intent(in) :: cell_idx
    integer(accs_int), intent(in) :: cell_face_ctr
  end subroutine set_face_location

  module subroutine set_neighbour_location(neighbour_location, cell_location, cell_neighbour_ctr)
    type(neighbour_locator), intent(out) :: neighbour_location
    type(cell_locator), intent(in) :: cell_location
    integer(accs_int), intent(in) :: cell_neighbour_ctr
  end subroutine set_neighbour_location
  
  end interface

  contains

  module subroutine set_global_matrix_size(mat, rows, columns, nnz, par_env)
    type(matrix_init_data), intent(inout) :: mat
    integer(accs_int), intent(in) :: rows
    integer(accs_int), intent(in) :: columns
    integer(accs_int), intent(in) :: nnz
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    mat%rglob = rows
    mat%cglob = columns
    mat%rloc = -1
    mat%cloc = -1
    mat%nnz = nnz
    mat%par_env => par_env
  end subroutine set_global_matrix_size

  module subroutine set_local_matrix_size(mat, rows, columns, nnz, par_env)
    type(matrix_init_data), intent(inout) :: mat
    integer(accs_int), intent(in) :: rows
    integer(accs_int), intent(in) :: columns
    integer(accs_int), intent(in) :: nnz
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    mat%rloc = rows
    mat%cloc = columns
    mat%rglob = -1
    mat%cglob = -1
    mat%nnz = nnz
    mat%par_env => par_env
  end subroutine set_local_matrix_size

  module subroutine set_face_location(face_location, geometry, cell_idx, cell_face_ctr)
    type(face_locator), intent(out) :: face_location
    type(mesh), target, intent(in) :: geometry
    integer(accs_int), intent(in) :: cell_idx
    integer(accs_int), intent(in) :: cell_face_ctr

    face_location%mesh => geometry
    face_location%cell_idx = cell_idx
    face_location%cell_face_ctr = cell_face_ctr
  end subroutine set_face_location

  module subroutine set_cell_location(cell_location, geometry, cell_idx)
    type(cell_locator), intent(out) :: cell_location
    type(mesh), target, intent(in) :: geometry
    integer(accs_int), intent(in) :: cell_idx

    ! XXX: Potentially expensive...
    if (cell_idx > size(geometry%idx_global)) then
      print *, "ERROR: trying to access cell I don't own!", cell_idx, geometry%nlocal
      stop
    else
      cell_location%mesh => geometry
      cell_location%cell_idx = cell_idx
    end if

  end subroutine set_cell_location

  module subroutine set_neighbour_location(neighbour_location, cell_location, cell_neighbour_ctr)
    type(neighbour_locator), intent(out) :: neighbour_location
    type(cell_locator), intent(in) :: cell_location
    integer(accs_int), intent(in) :: cell_neighbour_ctr

    ! integer(accs_int) :: nnb
    
    neighbour_location%mesh => cell_location%mesh
    neighbour_location%cell_idx = cell_location%cell_idx

    !! XXX: Safe, but would create a circular dependency...
    !! ! XXX: Potentially expensive...
    !! call count_neighbours(cell_location, nnb)
    !! if (cell_neighbour_ctr > nnb) then
    !!   print *, "ERROR: cell has fewer neighbours than neighbour count requested!"
    !!   stop
    !! else if (cell_neighbour_ctr < 1) then
    !!   print *, "ERROR: cell neighbour counter must be >= 1!"
    !! else
    !!   neighbour_location%cell_neighbour_ctr = cell_neighbour_ctr
    !! end if
    
    neighbour_location%cell_neighbour_ctr = cell_neighbour_ctr

    associate(mymesh => neighbour_location%mesh, &
         i => neighbour_location%cell_idx, &
         j => neighbour_location%cell_neighbour_ctr)
      if (mymesh%nbidx(j, i) == i) then
        print *, "ERROR: trying to set self as neighbour! Cell: ", i, j
      end if
    end associate
    
  end subroutine set_neighbour_location

end module types
