!>  Submodule file parallel_errors_mpi.smod
!
!>  @build mpi
!
!>  Implementation (using MPI) of parallel error handling

submodule (parallel) parallel_errors_mpi

  use mpi
  use parallel_types_mpi, only: parallel_environment_mpi

  implicit none

  contains

  !>  Error handling for parallel environment that uses MPI only
  module subroutine error_handling(error_code, error_category, par_env)

    integer, intent(in) :: error_code                   !< Variable the holds the error code
    character (len=*), intent (in) :: error_category    !< String description of the error category
    class(parallel_environment), intent(in) :: par_env  !< par_env
    
    integer(ccs_int) :: length
    integer :: ierr 
    character(len = MPI_MAX_ERROR_STRING) :: error_message 

    select type (par_env)

      type is (parallel_environment_mpi)

        if(error_category == "mpi") then
      
          if (error_code /= MPI_SUCCESS ) then
            call mpi_error_string(error_code, error_message, length, ierr)
            write(*,*) error_message
            call mpi_abort(par_env%comm, error_code, ierr)
          end if

        else 
          write(*,*) "Unsupported error categroy" 
          stop 1
        end if

      class default
        write(*,*) "Unsupported parallel environment"
        stop 1

    end select

  end subroutine

end submodule parallel_errors_mpi
