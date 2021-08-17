!> @brief Submodule file parallel_collectives_mpi.smod
!>
!> @details Implementation of parallel collectives using MPI

submodule (parallel) parallel_collectives_mpi

  use mpi

  implicit none

  contains

  !> @brief Global sum of integer scalars
  !>
  !> @param[in] integer input - Variable to be summed over on all ranks
  !> @param[out] integer result - Variable to hold the sum on all ranks
  !> @param[in] integer op - Variable that holds the reduction operation type
  !> @param[in] parallel_environment_mpi par_env
  module subroutine allreduce_scalar(input, result, op, par_env)

    class(*), intent(in) :: input
    class(*), intent(out) :: result
    integer, intent(in) :: op
    class(parallel_environment), intent(in) :: par_env
    integer :: ierr

    select type (par_env)
      type is (parallel_environment_mpi)   

      select type (input)
        type is (integer)
          call MPI_Allreduce(input, result, 1, MPI_INTEGER, op, par_env%comm, ierr)
        type is (double precision)
          call MPI_Allreduce(input, result, 1, MPI_DOUBLE, op, par_env%comm, ierr)
      class default
        write(*,*) "Unknown input data type"    
      end select

      call error_handling(ierr, par_env)
  
      class default
        write(*,*) "Unknown parallel environment"

    end select

  end subroutine

end submodule parallel_collectives_mpi