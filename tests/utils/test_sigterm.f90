!v A simple test to see if the basic functionality of stop_test is in place.
program test_sigterm
#include "ccs_macros.inc"
  use testing_lib
  use signal_handler, only: create_signal_handler, sigterm_issued

  implicit none
  integer :: i, j
  real :: a, b

  call init()

  sigterm_issued = .false.
  print *, "start", sigterm_issued
  call create_signal_handler()

  ! This is to avoid using a call to sleep, a subroutine from GNU extensions
  a = 0
  do i=0, 50000
    do j=0, 50000
      call random_number(b)
      a = a + b
      if (sigterm_issued) then
        exit
      end if
    end do
    if (sigterm_issued) then
      exit
    end if
  end do

  print *, "before check", sigterm_issued
  if (.not. sigterm_issued) then
    print *, "during check", sigterm_issued
    call error_abort("sigterm not called")
  end if
  print *, "after check", sigterm_issued

end program test_sigterm
