!v signal handler
module signal_handler
#include "ccs_macros.inc"

  use kinds
  use, intrinsic :: iso_c_binding, only : c_int, c_funloc, c_funptr

  implicit none

  private

  integer(ccs_int) :: signum
  logical, public :: sigterm_issued = .false.

  public :: create_signal_handler

  interface
    function c_signal(sig, func) bind(c, name='signal')
      import :: c_funptr, c_int
      implicit none
      integer(kind=c_int), intent(in), value :: sig
      type(c_funptr),      intent(in), value :: func
      type(c_funptr)                         :: c_signal
    end function c_signal
  end interface

contains

  subroutine create_signal_handler
    import :: c_funloc, c_funptr
    type(c_funptr) :: ptr_sigterm
    ! type(c_funptr) :: ptr_sigint
    integer(c_int), parameter :: SIGTERM = 15
    ! integer(c_int), parameter :: SIGINT = 2
    ptr_sigterm = c_signal(SIGTERM, c_funloc(catch_signal))
    ! ptr_siting = c_signal(SIGINT, c_funloc(catch_signal))
  end subroutine create_signal_handler

  subroutine catch_signal(signum) bind(c)
    integer(kind=c_int), intent(in), value :: signum
    sigterm_issued = .true.
  end subroutine catch_signal

end module signal_handler
