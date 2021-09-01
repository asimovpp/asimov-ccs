!> @brief Module file parallel_types.mod
!>
!> @details Module that defines the parallel environment types for ASiMoV-CCS
module parallel_types

  use mpi_f08

  implicit none

  private 

  !> @brief placeholder reduction operator type
  type, public :: reduction_operator
  end type reduction_operator

  !> @brief reduction operator type for MPI
  !>
  !> @details reduction operator type from MPI that holds
  !> the MPI operator values that are passed to reductions
  type, extends(reduction_operator), public :: reduction_operator_mpi
    type(mpi_op) :: sum_op
    type(mpi_op) :: min_op
    type(mpi_op) :: max_op
    type(mpi_op) :: prod_op
    type(mpi_op) :: land_op
    type(mpi_op) :: lor_op
    type(mpi_op) :: band_op
    type(mpi_op) :: bor_op
    type(mpi_op) :: maxloc_op
    type(mpi_op) :: minloc_op
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
    type(reduction_operator_mpi) :: rop
  end type parallel_environment_mpi

  !> @brief parallel environment type for CAF
  !>
  !> @details parallel environment type from CAF that holds
  !> image ID and the number of images
  type, extends(parallel_environment), public :: parallel_environment_caf
  end type parallel_environment_caf

  interface
    !> @brief Set the values of the reduction operators
    module subroutine set_reduction_operators(rop)
      class(reduction_operator), intent(inout) :: rop
    end subroutine
  end interface

  public :: set_reduction_operators

end module parallel_types