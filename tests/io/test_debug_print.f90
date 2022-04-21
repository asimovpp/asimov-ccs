program test_debug_print
#include "ccs_macros.inc"
  use utils, only : str
  use kinds, only : ccs_int, ccs_real
  implicit none

  call dprint("The numbers are " // str(24_ccs_int) // " and " // str(42_ccs_int, "(I3)") // " and " // str(6.28_ccs_real, "(F5.3)") // " and " // str(3.14_ccs_real)) 
end program test_debug_print
