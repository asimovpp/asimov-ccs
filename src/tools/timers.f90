!v Set of tools used to analyse discretisation errors
module timers
#include "ccs_macros.inc"

  use kinds
  use types
  use parallel_types, only: parallel_environment
  use parallel, only: timer

  use utils, only: exit_print

  implicit none

  private

  logical, dimension(:), allocatable :: has_started
  double precision, dimension(:), allocatable :: ticks
  double precision, dimension(:), allocatable :: tocks
  character(len=64), dimension(:), allocatable :: timer_names
  logical :: initialised = .false.

  public :: timer_init
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
      allocate(timer_names(0))
      initialised = .true.
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
  subroutine timer_register_start(timer_name, timer_index)
    character(len=*), intent(in) :: timer_name
    integer(ccs_int), intent(out) :: timer_index

    call timer_register(timer_name, timer_index)
    call timer_start(timer_index)

  end subroutine

  !> Register a new timer
  subroutine timer_register(timer_name, timer_index)
    character(len=*), intent(in) :: timer_name
    integer(ccs_int), intent(out) :: timer_index

    call timer_init()
    call timer_get_index(timer_name, timer_index)
    if (timer_index == -1) then
      timer_index = size(ticks) + 1
      timer_names = (/ timer_names, trim(timer_name) /)
      ticks = (/ ticks, 0.d0 /)
      tocks = (/ tocks, 0.d0 /)
      has_started = (/ has_started, .false. /)
    end if
    
  end subroutine

  !> Start a timer
  subroutine timer_start(timer_index)
    integer(ccs_int), intent(in) :: timer_index
    double precision :: tick
    double precision :: time

    call timer(tick)

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
    double precision :: time

    call timer_get_time(timer_index, time)

    if (par_env%proc_id == par_env%root) then
      ! repeat ' ' to right align
      write(*,'(A30, F10.4, A)')  repeat(' ', max(0, 29-len(trim(timer_names(timer_index))))) // trim(timer_names(timer_index)) // ":", time, " s"
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
