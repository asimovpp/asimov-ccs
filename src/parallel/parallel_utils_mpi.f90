!> @brief Submodule file parallel_utils_mpi.smod
!
!> @build mpi
!
!> @details Implementation (using MPI) of the parallel utilities

submodule (parallel) parallel_utils_mpi

  use mpi
  use parallel_types_mpi, only: parallel_environment_mpi

  implicit none

  contains

  !> @brief Synchronise the parallel environment
  !
  !> @param[in] parallel_environment_mpi par_env
  module subroutine sync(par_env)
  
    class(parallel_environment), intent(in) :: par_env
    integer :: ierr !> Error code

    select type (par_env)
      type is (parallel_environment_mpi)
         
        call mpi_barrier(par_env%comm, ierr)
        call error_handling(ierr, "mpi", par_env)

      class default
        write(*,*) "Unsupported parallel environment"

    end select

  end subroutine

  !> @brief Timer for MPI parallel environment
  !
  !> @param[inout] double precision tick - Variable that is assigned the current time
  module subroutine timer(tick)

    double precision, intent(out) :: tick

    tick = mpi_wtime()

  end subroutine

end submodule parallel_utils_mpi
