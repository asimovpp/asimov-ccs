submodule(profiler) profiler_caliper

  use caliper_mod

  implicit none

  type(ConfigManager) :: mgr

contains

  !> Initialise Caliper profiler
  module subroutine profiler_init()

    integer :: argc
    character(len=256) :: arg
    logical :: ret
    character(len=:), allocatable :: errmsg

    mgr = ConfigManager_new()
    do argc = 1, command_argument_count()
      if (argc .ge. 1) then
        call get_command_argument(argc, arg)
        if (arg(1:6) == '--cali') then
          call get_command_argument(argc + 1, arg)
          call mgr%add(arg)
          ret = mgr%error()
          if (ret) then
            errmsg = mgr%error_msg()
            write (*, *) 'ConfigManager: ', errmsg
          end if
        end if
      end if
    end do

    call mgr%start

  end subroutine

  !> Shutdown Caliper profiler
  module subroutine profiler_shutdown(par_env)
    class(parallel_environment), intent(in) :: par_env

    associate (foo => par_env)
    end associate
    call mgr%flush
    call cali_flush(0)
    call configmanager_delete(mgr)

  end subroutine

  !> Profile region -- start
  module subroutine profiler_begin_region(region_name)
    character(len=*), intent(in) :: region_name

    call cali_begin_region(region_name)

  end subroutine

  !> Profile region -- stop
  module subroutine profiler_end_region(region_name)
    character(len=*), intent(in) :: region_name

    call cali_end_region(region_name)

  end subroutine

end submodule
