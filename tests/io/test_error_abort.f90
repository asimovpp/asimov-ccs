program test_error_abort
#include "ccs_macros.inc"
  use testing_lib
  use utils, only: exit_print
  implicit none

  call init()
  call error_abort("Ending test")
  call fin()
end program test_error_abort
