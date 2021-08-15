!> @brief Submodule file parallel_utils_mpi.smod
!>
!> @details Implementation (using MPI) of the parallel utilities

submodule (parallel) parallel_utils_mpi

use mpi

implicit none

contains

  !> @brief Synchronise the parallel environment
  !>
  !> @param[in] iteger comm - MPI communicator to be synchronised
  module subroutine sync(comm)
  
    integer, intent(in) :: comm
    integer :: ierr, length, tmp_ierr
    character(len = MPI_MAX_ERROR_STRING) :: error_message

    call mpi_barrier(comm, ierr)
    call error_handling(ierr, comm)

  end subroutine

  !> @brief Timer for MPI parallel environment
  !>
  !> @param[inout] double precision tick - Variable that is assigned the current time
  module subroutine timer(tick)

    double precision, intent(inout) :: tick

    tick = mpi_wtime()

  end subroutine

  !> @brief Error handling for parallel environment
  !>
  !> @param[in] integer error_code - Variable the holds the error code
  !> @param[in] integer comm - Communicator to use when calling mpi_abort
  module subroutine error_handling(error_code, comm)

    integer, intent(in) :: error_code
    integer, intent(in) :: comm
    integer :: length, ierr
    character(len = MPI_MAX_ERROR_STRING) :: error_message

    if (error_code /= MPI_SUCCESS ) then
      call mpi_error_string(error_code, error_message, length, ierr)
      write(*,*) error_message(1:length)
      call mpi_abort(comm, error_code, ierr)
    end if

  end subroutine

end submodule parallel_utils_mpi