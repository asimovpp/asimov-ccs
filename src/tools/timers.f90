!v Set of tools used to analyse discretisation errors
module timers
#include "ccs_macros.inc"

  use kinds
  use types
  use parallel_types, only: parallel_environment
  use parallel, only: timer, is_root

  use utils, only: exit_print

  implicit none

  private

  logical, dimension(:), allocatable :: has_started
  double precision, dimension(:), allocatable :: ticks
  double precision, dimension(:), allocatable :: tocks
  integer(ccs_int), dimension(:), allocatable :: counters
  character(len=64), dimension(:), allocatable :: timer_names
  integer(ccs_int) :: total_index
  logical :: initialised = .false.

  public :: timer_init
  public :: timer_reset
  public :: timer_get_index
  public :: timer_register_start
  public :: timer_register
  public :: timer_start
  public :: timer_stop
  public :: timer_print
  public :: timer_print_all
  public :: timer_get_time

contains

  !> Initialise timer module
  subroutine timer_init()

    if (.not. initialised) then
      allocate(has_started(0))
      allocate(ticks(0))
      allocate(tocks(0))
      allocate(counters(0))
      allocate(timer_names(0))
      total_index = 0
      initialised = .true.
    end if

  end subroutine

  subroutine timer_reset()

    if (initialised) then
      deallocate(has_started)
      deallocate(ticks)
      deallocate(tocks)
      deallocate(counters)
      deallocate(timer_names)
      total_index = 0
      initialised = .false.
    end if

  end subroutine
  
  !> Get timer index from its name
  subroutine timer_get_index(timer_name, timer_index)
    character(len=*), intent(in) :: timer_name
    integer(ccs_int), intent(out) :: timer_index
    integer(ccs_int) :: i

    timer_index = -1

    do i=1, size(timer_names)
      if (trim(timer_name) == trim(timer_names(i))) then
        timer_index = i
        exit
      end if
    end do

  end subroutine

  !> Register and start a timer
  subroutine timer_register_start(timer_name, timer_index, is_total_time)
    character(len=*), intent(in) :: timer_name
    integer(ccs_int), intent(out) :: timer_index
    logical, optional, intent(in) :: is_total_time

    call timer_register(timer_name, timer_index, is_total_time=is_total_time)
    call timer_start(timer_index)

  end subroutine

  !> Register a new timer
  subroutine timer_register(timer_name, timer_index, is_total_time)
    character(len=*), intent(in) :: timer_name
    integer(ccs_int), intent(out) :: timer_index
    logical, optional, intent(in) :: is_total_time

    character(len=64) :: timer_name_local

    timer_name_local = timer_name

    call timer_init()
    call timer_get_index(timer_name_local, timer_index)
    if (timer_index == -1) then
      timer_index = size(ticks) + 1
      timer_names = [ timer_names, timer_name_local ]
      ticks = [ ticks, 0.d0 ]
      tocks = [ tocks, 0.d0 ]
      counters = [ counters, 0 ]
      has_started = [ has_started, .false. ]
      if (present(is_total_time)) then
        if (is_total_time) then
          total_index = timer_index
        end if
      end if
    end if
    
  end subroutine

  !> Start a timer
  subroutine timer_start(timer_index)
    integer(ccs_int), intent(in) :: timer_index
    double precision :: tick
    double precision :: time

    call timer(tick)
    counters(timer_index) = counters(timer_index) + 1

    ! Allows to start/stop a timer
    if (has_started(timer_index)) then
      call timer_get_time(timer_index, time)
      tick = tick - time
    else
      has_started(timer_index) = .true.
    end if

    ticks(timer_index) = tick

  end subroutine

  !> Stop a timer
  subroutine timer_stop(timer_index)
    integer(ccs_int), intent(in) :: timer_index

    call timer(tocks(timer_index))

  end subroutine

  !> Print time on rank 0
  subroutine timer_print(par_env, timer_index)
    class(parallel_environment), intent(in) :: par_env
    integer(ccs_int), intent(in) :: timer_index
    double precision :: time, total_time

    call timer_get_time(timer_index, time)

    if (is_root(par_env)) then
      ! repeat ' ' to right align
      write(*,'(A30, F10.4, A)', advance="no")  repeat(' ', max(0, 29-len(trim(timer_names(timer_index))))) // trim(timer_names(timer_index)) // ":", time, " s"

      if (total_index /= 0) then
        call timer_get_time(total_index, total_time)
        write(*, '(F10.2, A)', advance="no") 100*time / total_time, " %"
      end if

      if (counters(timer_index) >= 2) then
        write(*, '(F10.4, A, I10, A)', advance="no") time / counters(timer_index), " s/call", counters(timer_index), " calls"
      end if
      write(*, '(A)') ""
    end if

  end subroutine

  !> Print all the timers
  subroutine timer_print_all(par_env)
    class(parallel_environment), intent(in) :: par_env
    integer(ccs_int) :: timer_index

    do timer_index=1, size(ticks)
      call timer_print(par_env, timer_index)
    end do

  end subroutine

  !> Get raw time from timer
  subroutine timer_get_time(timer_index, time)
    integer(ccs_int), intent(in) :: timer_index
    double precision :: time

    time = tocks(timer_index) - ticks(timer_index)

  end subroutine

end module timers
