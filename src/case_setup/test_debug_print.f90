program test_debug_print
#include "ccs_macros.inc"
  use utils, only : str
  use kinds, only : ccs_int, ccs_real
  implicit none

  call dprint("The numbers are " // trim(adjustl(str(42_ccs_int))) // " and also " // trim(adjustl(str(3.14_ccs_real, '(F5.3)')))) 
  call dprint("Other numbers are " // str(24_ccs_int) // " and also " // str(6.28_ccs_real)) 
end program test_debug_print
