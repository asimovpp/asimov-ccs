#include "ccs_macros.h"
program test_debug_print
  use utils, only : debug_print, str
  use kinds, only : ccs_int, ccs_real
  implicit none

  call dprint("The numbers are " // trim(adjustl(str(42_ccs_int))) // " and also " // trim(adjustl(str(3.14_ccs_real, '(F5.3)')))) 
end program test_debug_print
