!> @brief Module file parallel.mod
!>
!> @details Module that defines the parallel interace for ASiMoV-CCS

module parallel

  use parallel_types

  implicit none

  private

  interface

    !> @brief Create the parallel environment
    module subroutine setup_parallel_environment(par_env)
      class(parallel_environment), intent(out) :: par_env
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
    module subroutine allreduce_scalar(input, result, op, par_env)
      class(*), intent(in) :: input
      class(*), intent(out) :: result
      integer, intent(in) :: op
      class(parallel_environment), intent(in) :: par_env
    end subroutine

    !> @brief Error handling for parallel environment
    module subroutine error_handling(error_code, par_env)
      integer, intent(in) :: error_code
      class(parallel_environment), intent(in) :: par_env
    end subroutine

  end interface

  interface allreduce
    module procedure allreduce_scalar
  end interface allreduce

  public :: setup_parallel_environment
  public :: cleanup_parallel_environment
  public :: sync
  public :: timer
  public :: allreduce

end module parallel