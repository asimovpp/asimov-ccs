!v Profiler module
module profiler

  use parallel_types, only: parallel_environment

  implicit none

  private

  public :: profiler_init
  public :: profiler_shutdown
  public :: profiler_begin_region
  public :: profiler_end_region

  interface

    !> Initialise profiler
    module subroutine profiler_init()
    end subroutine

    !> Shutdown profiler
    module subroutine profiler_shutdown(par_env)
      class(parallel_environment), intent(in) :: par_env
    end subroutine

    !> Profile region -- start
    module subroutine profiler_begin_region(region_name)
      character(len=*), intent(in) :: region_name
    end subroutine

    !> Profile region -- stop
    module subroutine profiler_end_region(region_name)
      character(len=*), intent(in) :: region_name
    end subroutine

  end interface

end module profiler
