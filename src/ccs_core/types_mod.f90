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
    real(accs_real) :: h, Af
    real(accs_real), dimension(:), allocatable :: vol        !> Cell volumes
    real(accs_real), dimension(:, :), allocatable :: xc      !> Cell centres (dimension, cell)
    real(accs_real), dimension(:, :, :), allocatable :: xf   !> Face centres (dimension, face, cell)
  end type mesh

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

end module types
