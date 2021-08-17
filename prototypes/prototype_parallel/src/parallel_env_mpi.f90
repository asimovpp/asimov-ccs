!> @brief Submodule file parallel_env_mpi.smod
!>
!> @details Implementation of the parallel environment using MPI

submodule (parallel) parallel_env_mpi

  use mpi

  implicit none

  contains

  !> @brief Create the MPI parallel environment
  !>
  !> @param[out] parallel_environment_mpi par_env
  module subroutine setup_parallel_environment(par_env)

    integer :: ierr, length, tmp_ierr
    character(len = MPI_MAX_ERROR_STRING) :: error_message

    class(parallel_environment), intent(out) :: par_env

    select type (par_env)
      type is (parallel_environment_mpi)   
        par_env%comm = MPI_COMM_WORLD
        call mpi_init(ierr)
        call error_handling(ierr, par_env)

        call mpi_comm_rank(par_env%comm, par_env%rank, ierr)
        call error_handling(ierr, par_env)

        call mpi_comm_size(par_env%comm, par_env%numprocs, ierr)
        call error_handling(ierr, par_env)
   end select

  end subroutine

  !> @brief Cleanup the MPI parallel environment
  !>
  !> @param[in] parallel_environment_mpi par_env
  module subroutine cleanup_parallel_environment(par_env)

    class(parallel_environment), intent(in) :: par_env
    integer :: ierr, length, tmp_ierr
    character(len = MPI_MAX_ERROR_STRING) :: error_message

    select type (par_env)
      type is (parallel_environment_mpi)   
      call mpi_finalize(ierr)
      call error_handling(ierr, par_env)
    end select

    end subroutine

  end submodule parallel_env_mpi