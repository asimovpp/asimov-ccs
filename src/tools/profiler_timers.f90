!v Interface to timer-based profiler
submodule(profiler) profiler_timers

  implicit none

contains

  !> Initialise timer-based profiler
  module subroutine profiler_init()
    use timers, only: timer_init

    call timer_init()

  end subroutine

  !> Shutdown timer-based profiler
  module subroutine profiler_shutdown(par_env)
    use timers, only: timer_print_all, timer_export_csv, timer_reset
    use parallel_types_mpi, only: parallel_environment_mpi
  
    class(parallel_environment), intent(in) :: par_env
    
    select type(par_env)
      type is (parallel_environment_mpi)
        call timer_export_csv(par_env)
        call timer_print_all(par_env)
      class default
        error stop "Unknown parallel environment"
    end select

    call timer_reset()

  end subroutine

  module subroutine profiler_begin_region(region_name)
    use timers, only: timer_register_start

    character(len=*), intent(in) :: region_name
    integer :: timer_index

    call timer_register_start(region_name, timer_index)

  end subroutine

  module subroutine profiler_end_region(region_name)
    use timers, only: timer_stop, timer_get_index

    character(len=*), intent(in) :: region_name
    integer :: timer_index

    call timer_get_index(region_name, timer_index)
    call timer_stop(timer_index)

  end subroutine


end submodule
