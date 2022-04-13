#include "ccs_macros.h"
program test_debug_print
  use utils, only : diag_null, diag_print, str
  use kinds, only : ccs_int, ccs_real
  implicit none

  call diag("The numbers are " // str(42_ccs_int) // " and also " // str(3.14_ccs_real)) 
end program test_debug_print
