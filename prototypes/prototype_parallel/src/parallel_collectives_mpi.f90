!> @brief Submodule file parallel_collectives_mpi.smod
!>
!> @details Implementation of parallel collectives using MPI

submodule (parallel) parallel_collectives_mpi

  use mpi_f08

  implicit none

  contains

  !> @brief Global sum of integer scalars
  !>
  !> @param[in] integer input_value - Variable to be summed over on all ranks
  !> @param[out] integer result_value - Variable to hold the sum on all ranks
  !> @param[in] integer op - Variable that holds the reduction operation type
  !> @param[in] parallel_environment_mpi par_env
  module subroutine allreduce_scalar(input_value, result_value, rop, par_env)

    class(*), intent(in) :: input_value
    class(*), intent(out) :: result_value
    class(reduction_operator), intent(in) :: rop
    class(parallel_environment), intent(in) :: par_env
    integer :: ierr

    select type (par_env)

    type is (parallel_environment_mpi)   

      select type(rop)
        type is(reduction_operator_mpi)

!        write(*,*) "Type is correct"
        select type (input_value)
        type is (integer)
          call MPI_Allreduce(input_value, result_value, 1, MPI_INTEGER, rop%op, par_env%comm, ierr)

        type is (double precision)
          call MPI_Allreduce(input_value, result_value, 1, MPI_DOUBLE, rop%op, par_env%comm, ierr)

        class default
          write(*,*) "Unsupported input data type"    

      end select

      class default
        write(*,*) "Unsupported reduction operation type"

      end select

      call error_handling(ierr, par_env)
  
    class default
      write(*,*) "Unsupported parallel environment"

    end select

  end subroutine

end submodule parallel_collectives_mpi