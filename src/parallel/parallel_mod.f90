!> @brief Module file parallel.mod
!
!> @details Module that defines the parallel interace for ASiMoV-CCS

module parallel

  use parallel_types
  use mpi_f08

  implicit none

  private

  interface

    !> @brief Create the parallel environment
    module subroutine initialise_parallel_environment(par_env)
      class(parallel_environment), allocatable, intent(out) :: par_env
    end subroutine

    !> @brief Cleanup the parallel environment
    module subroutine cleanup_parallel_environment(par_env)
      class(parallel_environment), intent(in) :: par_env
    end subroutine

    !> @brief Synchronise the parallel environment
    module subroutine sync(par_env)
      class(parallel_environment), intent(in) :: par_env
    end subroutine

    !> @brief Timer for parallel environment
    module subroutine timer(tick)
      double precision, intent(out) :: tick
    end subroutine

    !> @brief Global reduction of integer scalars
    module subroutine allreduce_scalar(input_value, result_value, rop, par_env)
      class(*), intent(in) :: input_value
      class(*), intent(inout) :: result_value
      class(reduction_operator), intent(in) :: rop
      class(parallel_environment), intent(in) :: par_env
    end subroutine

    !> @brief Error handling for parallel environment
    module subroutine error_handling(error_code, error_cat, par_env)
      integer, intent(in) :: error_code
      character (len=*), intent (in) :: error_cat
      class(parallel_environment), intent(in) :: par_env
    end subroutine

  end interface

  interface allreduce
    module procedure allreduce_scalar
  end interface allreduce

  public :: initialise_parallel_environment
  public :: cleanup_parallel_environment
  public :: sync
  public :: timer
  public :: allreduce

end module parallel
