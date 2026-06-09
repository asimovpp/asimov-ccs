!v signal handler
module signal_handler
#include "ccs_macros.inc"

  use kinds
  use, intrinsic :: iso_c_binding, only : c_int, c_funloc, c_funptr

  implicit none

  private

  logical, public :: sigterm_issued = .false.

  public :: create_signal_handler

  !> Fortran interface for the C function `signal`
  interface
    function c_signal(sig, func) bind(c, name='signal')
      import :: c_funptr, c_int
      implicit none
      integer(kind=c_int), intent(in), value :: sig      !< the POSIX signal number to be captured
      type(c_funptr),      intent(in), value :: func     !< the signal handler function to be called
                                                         !< when the target signal is captured
      type(c_funptr)                         :: c_signal !< this function
    end function c_signal
  end interface

contains

  !> Create signal handler to pass to c_signal
  subroutine create_signal_handler
    type(c_funptr) :: ptr_sigterm             !< function pointer to c_signal function that creates the signal trap
    integer(c_int), parameter :: SIGTERM = 15 !< POSIX signal number to trap
    ptr_sigterm = c_signal(SIGTERM, c_funloc(catch_signal))
  end subroutine create_signal_handler

  !> Routine that gets called when signal `signum` is trapped
  subroutine catch_signal(signum) bind(c)
    integer(kind=c_int), intent(in), value :: signum !< POSIX signal number that, once captured, triggers this function
    associate(foo => signum)
    end associate
    sigterm_issued = .true. !< this public argument is checked by the core_solver and, when true, triggers a stop
  end subroutine catch_signal

end module signal_handler
