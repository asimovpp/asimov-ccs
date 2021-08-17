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
  module subroutine allreduce_integer(input, result, op, par_env)

    integer, intent(in) :: input
    integer, intent(out) :: result
    integer, intent(in) :: op
    class(parallel_environment), intent(in) :: par_env
    integer :: ierr

    select type (par_env)
      type is (parallel_environment_mpi)   

        call MPI_Allreduce(input, result, 1, MPI_INTEGER, op, par_env%comm, ierr)
        call error_handling(ierr, par_env)

    end select

  end subroutine

  !> @brief Global sum of double precision real scalars
  !>
  !> @param[in] double precision input - Variable to be summed over on all ranks
  !> @param[out] double precision result - Variable to hold the sum on all ranks
  !> @param[in] integer op - Variable that holds the reduction operation type
  !> @param[in] parallel_environment_mpi par_env
  module subroutine allreduce_double(input, result, op, par_env)

    double precision, intent(in) :: input
    double precision, intent(out) :: result
    integer, intent(in) :: op
    class(parallel_environment), intent(in) :: par_env
    integer :: ierr

    select type (par_env)
      type is (parallel_environment_mpi)   

        call MPI_Allreduce(input, result, 1, MPI_DOUBLE, op, par_env%comm, ierr)
        call error_handling(ierr, par_env)

    end select
  end subroutine

end submodule parallel_collectives_mpi