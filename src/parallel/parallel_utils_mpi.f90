!> @brief Submodule file parallel_utils_mpi.smod
!>
!> @details Implementation (using MPI) of the parallel utilities

submodule (parallel) parallel_utils_mpi

  use mpi_f08
  use parallel_types_mpi, only: parallel_environment_mpi

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
      call error_handling(ierr, "mpi", par_env)

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
  !> @param[in] string error_cat - String description of the error category
  !> @param[in] parallel_environment_mpi par_env
  module subroutine error_handling(error_code, error_cat, par_env)

    integer, intent(in) :: error_code
    character (len=*), intent (in) :: error_cat
    class(parallel_environment), intent(in) :: par_env
    integer :: length, ierr
    character(len = MPI_MAX_ERROR_STRING) :: error_message

    select type (par_env)

    type is (parallel_environment_mpi)

    if(error_cat == "mpi") then
      
      if (error_code /= MPI_SUCCESS ) then
        call mpi_error_string(error_code, error_message, length, ierr)
        write(*,*) error_message(1:length)
        call mpi_abort(par_env%comm, error_code, ierr)
      end if

    else if (error_cat == "petsc") then

      if (error_code /= 0) then
        write(*,*) "PETSc error"
        call mpi_abort(par_env%comm, error_code, ierr)
      end if

    else 
      write(*,*) "Unsupported error categroy" 
    end if

    class default
      write(*,*) "Unsupported parallel environment"

    end select

  end subroutine

end submodule parallel_utils_mpi