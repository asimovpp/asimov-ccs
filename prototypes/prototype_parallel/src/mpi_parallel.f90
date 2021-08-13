!> @brief Submodule file mpi_parallel.smod
!>
!> @details Implementation (using MPI) of the parallel module interface

submodule (parallel) mpi_parallel

  use mpi

  implicit none

contains

  !> @brief Create the MPI parallel environment
  !>
  !> @param[out] integer comm - MPI communicator
  !> @param[out] integer rank - MPI rank
  !> @param[out] integer size - Total number of MPI ranks
  !> @param[out] integer ierr - error code
  module subroutine setup_parallel_environment(comm, rank, size, ierr)
    integer, intent(out) :: comm
    integer, intent(out) :: rank
    integer, intent(out) :: size
    integer, intent(out) :: ierr

    comm = MPI_COMM_WORLD
    call mpi_init(ierr)
    call mpi_comm_rank(comm, rank, ierr)
    call mpi_comm_size(comm, size, ierr)

  end subroutine

  !> @brief Cleanup the MPI parallel environment
  !>
  !> @param[out] integer ierr - error code
  module subroutine cleanup_parallel_environment(ierr)
    integer, intent(out) :: ierr

    call mpi_finalize(ierr)

  end subroutine

  end submodule mpi_parallel