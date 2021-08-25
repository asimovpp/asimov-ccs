!> @brief Module file parallel_types.mod
!>
!> @details Module that defines the parallel environment types for ASiMoV-CCS
module parallel_types

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
    integer(kind=4) :: sum_op
    integer(kind=4) :: min_op
    integer(kind=4) :: max_op
    integer(kind=4) :: prod_op
    integer(kind=4) :: land_op
    integer(kind=4) :: lor_op
    integer(kind=4) :: band_op
    integer(kind=4) :: bor_op
    integer(kind=4) :: maxloc_op
    integer(kind=4) :: minloc_op
  end type reduction_operator_mpi

  interface

    !> @brief Set the values of the reduction operators
    module subroutine set_reduction_operators(rop)
      class(reduction_operator), intent(out) :: rop
    end subroutine

  end interface

  !> @brief placeholder parallel environment type
  type, public :: parallel_environment
  end type parallel_environment

  !> @brief parallel environment type for MPI
  !>
  !> @details parallel environment type from MPI that holds
  !> rank, number of processes and communicator
  type, extends(parallel_environment), public :: parallel_environment_mpi
    integer :: rank
    integer :: numprocs
    integer :: comm
    type(reduction_operator_mpi) :: rop
  end type parallel_environment_mpi

  public :: set_reduction_operators

end module parallel_types