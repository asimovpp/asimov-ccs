!v Module file parallel_types.mod
!
!  Module that defines the parallel environment types for ASiMoV-CCS
!
!  @build mpi
module parallel_types_mpi

  use mpi
  use parallel_types, only: parallel_environment, reduction_operator

  implicit none

  private

  !v reduction operator type for MPI
  !
  !  reduction operator type from MPI that holds
  !  the MPI operator values that are passed to
  !  reductions
  type, extends(reduction_operator), public :: reduction_operator_mpi
    integer :: op
  end type reduction_operator_mpi

  !v parallel environment type for MPI
  !
  !  parallel environment type from MPI that holds
  !  a communicator and reduction operators in
  !  addition to the common parameters
  type, extends(parallel_environment), public :: parallel_environment_mpi
    integer :: comm
    integer :: sum_op
    integer :: min_op
    integer :: max_op
    integer :: prod_op
    integer :: land_op
    integer :: lor_op
    integer :: band_op
    integer :: bor_op
    integer :: maxloc_op
    integer :: minloc_op
  contains
    procedure :: set_rop => set_mpi_reduction_operator
  end type parallel_environment_mpi

contains
  !> Set the values of the reduction operators
  subroutine set_mpi_reduction_operator(this)
    class(parallel_environment_mpi), intent(inout) :: this
    this%sum_op = MPI_SUM
    this%min_op = MPI_MIN
    this%max_op = MPI_MAX
    this%prod_op = MPI_PROD
    this%land_op = MPI_LAND
    this%lor_op = MPI_LOR
    this%band_op = MPI_BAND
    this%bor_op = MPI_BOR
    this%maxloc_op = MPI_MAXLOC
    this%minloc_op = MPI_MINLOC
  end subroutine

end module parallel_types_mpi
