!v Interface to Caliper profiler
module profiler
#include "ccs_macros.inc"

  use caliper_mod
  use iso_c_binding, ONLY : C_INT64_T

  implicit none

  private

  type(ConfigManager) :: mgr

  public :: profiler_init
  public :: profiler_shutdown

contains

  !> Initialise Caliper profiler
  subroutine profiler_init()

    integer :: argc
    character(len=256)    :: arg
    logical               :: ret
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
            write(*,*) 'ConfigManager: ', errmsg
          endif
        endif
      endif
    end do

    call mgr%start

  end subroutine

  !> Shutdown Caliper profiler
  subroutine profiler_shutdown()

    call mgr%flush
    call configmanager_delete(mgr)

    end subroutine

end module profiler
