!> @brief Submodule file parallel_utils_caf.smod
!>
!> @details Implementation (using CAF) of the parallel utilities

submodule (parallel) parallel_utils_caf

  implicit none

  contains

  !> @brief Synchronise the parallel environment
  !>
  !> @param[in] parallel_environment_mpi par_env
  module subroutine sync(par_env)
  
    class(parallel_environment), intent(in) :: par_env

    select type (par_env)

    type is (parallel_environment_caf)   
      sync all

    class default
      write(*,*) "Unsupported parallel environment"

    end select

  end subroutine

  !> @brief Timer for CAF parallel environment
  !>
  !> @param[inout] double precision tick - Variable that is assigned the current time
  module subroutine timer(tick)

    double precision, intent(out) :: tick

    call cpu_time(tick)

  end subroutine

end submodule parallel_utils_caf