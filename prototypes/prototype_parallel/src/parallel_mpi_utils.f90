!> @brief Submodule file parallel_mpi_utils.smod
!>
!> @details Implementation (using MPI) of the parallel utils interface

submodule (parallel) parallel_mpi_utils

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

    if (ierr /= MPI_SUCCESS ) then
      call mpi_error_string(ierr, error_message, length, tmp_ierr)
      write(*,*) error_message(1:length)
      call mpi_abort(comm, ierr, tmp_ierr)
    end if

  end subroutine

  !> @brief Timer for MPI parallel environment
  !>
  !> @param[inout] double precision tick - Variable that is assigned the current time
  module subroutine timer(tick)

    double precision, intent(inout) :: tick

    tick = mpi_wtime()

  end subroutine

end submodule parallel_mpi_utils