!> @brief Module file parallel_types.mod
!>
!> @details Module that defines the parallel environment types for ASiMoV-CCS
module parallel_types

  use mpi_f08

  implicit none

  private 

  public :: set_mpi_reduction_operator

  !> @brief placeholder reduction operator type
  type, public :: reduction_operator
  end type reduction_operator

  !> @brief reduction operator type for MPI
  !>
  !> @details reduction operator type from MPI that holds
  !> the MPI operator values that are passed to reductions
  type, extends(reduction_operator), public :: reduction_operator_mpi
    type(mpi_op) :: op
  end type reduction_operator_mpi

  !> @brief parallel environment type with common parameters
  !> process id, number of processes and root process
  type, public :: parallel_environment
    integer :: proc_id
    integer :: num_procs
    integer :: root
  end type parallel_environment

  !> @brief parallel environment type for MPI
  !>
  !> @details parallel environment type from MPI that holds
  !> a communicator and reduction operators in addition
  !> to the common parameters
  type, extends(parallel_environment), public :: parallel_environment_mpi
    type(mpi_comm) :: comm
    type(reduction_operator_mpi) :: sum_op
    type(reduction_operator_mpi) :: min_op
    type(reduction_operator_mpi) :: max_op
    type(reduction_operator_mpi) :: prod_op
    type(reduction_operator_mpi) :: land_op
    type(reduction_operator_mpi) :: lor_op
    type(reduction_operator_mpi) :: band_op
    type(reduction_operator_mpi) :: bor_op
    type(reduction_operator_mpi) :: maxloc_op
    type(reduction_operator_mpi) :: minloc_op
    contains
      procedure :: set_rop => set_mpi_reduction_operator
  end type parallel_environment_mpi

  !> @brief parallel environment type for CAF
  !>
  !> @details parallel environment type from CAF that holds
  !> image ID and the number of images
  type, extends(parallel_environment), public :: parallel_environment_caf
  end type parallel_environment_caf

  interface
    !> @brief Set the values of the reduction operators
    module subroutine set_mpi_reduction_operator(this)
      class(parallel_environment), intent(inout) :: this
    end subroutine
  end interface


end module parallel_types