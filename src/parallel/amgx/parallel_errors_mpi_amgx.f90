!v Submodule file parallel_errors_mpi_petsc.smod
!
!  @details Implementation (using MPI and PETSc) of parallel error handling
!
!  @build mpi petsc

submodule(parallel) parallel_errors_mpi_petsc
#include "ccs_macros.inc"

  use utils, only: exit_print
  use mpi
  use parallel_types_mpi, only: parallel_environment_mpi

  implicit none

contains

  !v Error handling for parallel environment that use
  !  both MPI and PETSc
  module subroutine error_handling(error_code, error_category, par_env)

    integer, intent(in) :: error_code !< integer error_code - Variable the holds the error code
    character(len=*), intent(in) :: error_category !< string error_category - String description of the error category
    class(parallel_environment), intent(in) :: par_env !< parallel_environment_mpi

    integer(ccs_int) :: length
    integer :: ierr
    character(len=MPI_MAX_ERROR_STRING) :: error_message

    select type (par_env)

    type is (parallel_environment_mpi)

      if (error_category == "mpi") then

        if (error_code /= MPI_SUCCESS) then
          call mpi_error_string(error_code, error_message, length, ierr)
          write (*, *) error_message
          call mpi_abort(par_env%comm, error_code, ierr)
        end if

      else if (error_category == "petsc") then

        if (error_code /= 0) then
          write (*, *) "PETSc error: ", error_code
          call mpi_abort(par_env%comm, error_code, ierr)
        end if

      else
        call error_abort("Unsupported error categroy")
      end if

    class default
      call error_abort("Unsupported parallel environment")

    end select

  end subroutine

end submodule parallel_errors_mpi_petsc
