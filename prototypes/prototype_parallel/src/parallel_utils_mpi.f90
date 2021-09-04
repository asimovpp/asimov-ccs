!> @brief Submodule file parallel_utils_mpi.smod
!>
!> @details Implementation (using MPI) of the parallel utilities

submodule (parallel) parallel_utils_mpi

  use mpi_f08

  implicit none

  contains

  !> @brief Synchronise the parallel environment
  !>
  !> @param[in] parallel_environment_mpi par_env
  module subroutine sync(par_env)
  
    class(parallel_environment), intent(in) :: par_env
    integer :: ierr

    select type (par_env)

    type is (parallel_environment_mpi)   
      call mpi_barrier(par_env%comm, ierr)
      call error_handling(ierr, par_env)

    class default
      write(*,*) "Unsupported parallel environment"

    end select

  end subroutine

  !> @brief Timer for MPI parallel environment
  !>
  !> @param[inout] double precision tick - Variable that is assigned the current time
  module subroutine timer(tick)

    double precision, intent(out) :: tick

    tick = mpi_wtime()

  end subroutine

  !> @brief Error handling for parallel environment
  !>
  !> @param[in] integer error_code - Variable the holds the error code
  !> @param[in] parallel_environment_mpi par_env
  module subroutine error_handling(error_code, par_env)

    integer, intent(in) :: error_code
    class(parallel_environment), intent(in) :: par_env
    integer :: length, ierr
    character(len = MPI_MAX_ERROR_STRING) :: error_message

    select type (par_env)

    type is (parallel_environment_mpi)   
      if (error_code /= MPI_SUCCESS ) then
        call mpi_error_string(error_code, error_message, length, ierr)
        write(*,*) error_message(1:length)
        call mpi_abort(par_env%comm, error_code, ierr)
      end if

    class default
      write(*,*) "Unsupported parallel environment"

    end select

  end subroutine

end submodule parallel_utils_mpi