!> @brief Submodule file parallel_env_mpi.smod
!>
!> @details Implementation of the parallel environment using MPI

submodule (parallel) parallel_env_mpi

  use mpi

  implicit none

contains

  !> @brief Create the MPI parallel environment
  !>
  !> @param[out] integer comm - MPI communicator
  !> @param[out] integer rank - MPI rank
  !> @param[out] integer numprocs - Total number of MPI ranks
  module subroutine setup_parallel_environment(comm, rank, numprocs)

    integer, intent(out) :: comm
    integer, intent(out) :: rank
    integer, intent(out) :: numprocs
    integer :: ierr, length, tmp_ierr
    character(len = MPI_MAX_ERROR_STRING) :: error_message

    comm = MPI_COMM_WORLD
    call mpi_init(ierr)
    call error_handling(ierr, comm)

    call mpi_comm_rank(comm, rank, ierr)
    call error_handling(ierr, comm)

    call mpi_comm_size(comm, numprocs, ierr)
    call error_handling(ierr, comm)

  end subroutine

  !> @brief Cleanup the MPI parallel environment
  !>
  !> @param[in] integer comm - MPI communicator to be cleaned up
  module subroutine cleanup_parallel_environment(comm)

    integer, intent(in) :: comm
    integer :: ierr, length, tmp_ierr
    character(len = MPI_MAX_ERROR_STRING) :: error_message

    call mpi_finalize(ierr)
    call error_handling(ierr, comm)

    end subroutine

  end submodule parallel_env_mpi