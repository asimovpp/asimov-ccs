!> @brief Module file parallel_types.mod
!
!> @build mpi
!
!> @details Module that defines the parallel environment types for ASiMoV-CCS
module parallel_types_mpi

  use mpi_f08
  use parallel_types, only: parallel_environment, & 
                            reduction_operator

  implicit none

  private 

  !> @brief reduction operator type for MPI
  !
  !> @details reduction operator type from MPI that holds
  !!          the MPI operator values that are passed to
  !!          reductions
  type, extends(reduction_operator), public :: reduction_operator_mpi
    type(mpi_op) :: op
  end type reduction_operator_mpi

  !> @brief parallel environment type for MPI
  !
  !> @details parallel environment type from MPI that holds
  !!          a communicator and reduction operators in
  !!          addition to the common parameters
  type, extends(parallel_environment), public :: parallel_environment_mpi
    type(mpi_comm) :: comm
    integer :: comm_int
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

  contains 
  !> @brief Set the values of the reduction operators
  subroutine set_mpi_reduction_operator(this)
    class(parallel_environment_mpi), intent(inout) :: this
    this%sum_op%op = MPI_SUM
    this%min_op%op = MPI_MIN
    this%max_op%op = MPI_MAX
    this%prod_op%op = MPI_PROD
    this%land_op%op = MPI_LAND
    this%lor_op%op = MPI_LOR
    this%band_op%op = MPI_BAND
    this%bor_op%op = MPI_BOR
    this%maxloc_op%op = MPI_MAXLOC
    this%minloc_op%op = MPI_MINLOC
  end subroutine

end module parallel_types_mpi
