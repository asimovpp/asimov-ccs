!!!        FILE: mod_accs_test_utils.f90
!!! DESCRIPTION: Module providing utility functions/subroutines for the ASiMoV-CCS testing
!!!              framework.

module accs_test_utils

  use constants
  
  implicit none

  private

  public :: accs_test_scale

contains

  pure function accs_test_scale(input)

    real(accs_real), intent(in) :: input
    real(accs_real) :: accs_test_scale
    
    accs_test_scale = TESTSCALE * input

  end function accs_test_scale

end module accs_test_utils
