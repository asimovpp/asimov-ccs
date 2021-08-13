!> @brief Submodule file parallel_mpi.smod
!>
!> @details Implementation (using MPI) of the parallel module interface

submodule (parallel) parallel_mpi

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
    
    if (ierr /= MPI_SUCCESS ) then
      call mpi_error_string(ierr, error_message, length, tmp_ierr)
      write(*,*) error_message(1:length)
    end if

    call mpi_comm_rank(comm, rank, ierr)
    call mpi_comm_size(comm, numprocs, ierr)

  end subroutine

  !> @brief Cleanup the MPI parallel environment
  module subroutine cleanup_parallel_environment()

    integer :: ierr, length, tmp_ierr
    character(len = MPI_MAX_ERROR_STRING) :: error_message

    call mpi_finalize(ierr)

    if (ierr /= MPI_SUCCESS ) then
        call mpi_error_string(ierr, error_message, length, tmp_ierr)
        write(*,*) error_message(1:length)
    end if

    end subroutine

  end submodule parallel_mpi