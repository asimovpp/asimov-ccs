!> @brief Submodule file parallel_collectives_mpi.smod
!>
!> @details Implementation of parallel collectives using MPI

submodule (parallel) parallel_collectives_mpi

use mpi

implicit none

contains

  !> @brief Global sum of integer scalars
  !>
  !> @param[in] integer x - Variable to be summed over on all ranks
  !> @param[out] integer sum - Variable to hold the sum on all ranks
  !> @param[comm] integer comm - Communicator to perform the sum over
  module subroutine global_sum_integer(x, sum, comm)

    integer, intent(in) :: x
    integer, intent(out) :: sum
    integer, intent(in) :: comm
    integer :: ierr

    call MPI_Allreduce(x, sum, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
    call error_handling(ierr, comm)

  end subroutine

  !> @brief Global sum of double precision real scalars
  !>
  !> @param[in] double precision x - Variable to be summed over on all ranks
  !> @param[out] double precision sum - Variable to hold the sum on all ranks
  !> @param[comm] double precision comm - Communicator to perform the sum over
  module subroutine global_sum_double(x, sum, comm)

    double precision, intent(in) :: x
    double precision, intent(out) :: sum
    integer, intent(in) :: comm
    integer :: ierr

    call MPI_Allreduce(x, sum, 1, MPI_DOUBLE, MPI_SUM, comm, ierr)
    call error_handling(ierr, comm)
    
  end subroutine

end submodule parallel_collectives_mpi