!v Interface to Likwid-based profiler
submodule(profiler) profiler_likwid

  use likwid

  implicit none

contains

  !> Initialise timer-based profiler
  module subroutine profiler_init()
    call likwid_markerInit()
  end subroutine

  !> Shutdown timer-based profiler
  module subroutine profiler_shutdown(par_env)
    class(parallel_environment), intent(in) :: par_env

    associate (foo => par_env)
    end associate
    call likwid_markerClose()

  end subroutine

  module subroutine profiler_begin_region(region_name)
    character(len=*), intent(in) :: region_name

    call likwid_markerRegisterRegion(region_name)
    call likwid_markerStartRegion(region_name)

  end subroutine

  module subroutine profiler_end_region(region_name)
    character(len=*), intent(in) :: region_name

    call likwid_markerStopRegion(region_name)
  end subroutine

end submodule
