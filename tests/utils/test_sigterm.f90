!v A simple test to see if the basic functionality of stop_test is in place.
program test_sigterm
#include "ccs_macros.inc"
  use testing_lib
  use signal_handler, only: create_signal_handler, sigterm_issued
  use parallel, only: is_root
  use mpi

  implicit none
  integer :: i, j
  real :: a, b

  call init()

  call create_signal_handler()

  ! This is to avoid using a call to sleep, a subroutine from GNU extensions
  a = 0
  do i=0, 50000
    do j=0, 50000
      call random_number(b)
      a = a + b
    end do
    if (sigterm_issued) then
      exit
    end if
  end do

  if (is_root(par_env) .and. (.not. sigterm_issued)) then
    call error_abort("SIGTERM not caught")
  end if

  call fin()

end program test_sigterm
