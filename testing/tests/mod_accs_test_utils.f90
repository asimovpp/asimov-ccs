!!!        FILE: mod_accs_test_utils.f90
!!! DESCRIPTION: Module providing utility functions/subroutines for the ASiMoV-CCS testing
!!!              framework.

module accs_test_utils

  use constants, only : accs_real, accs_real_eps
  
  implicit none

  private

  public :: accs_test_scale, accs_test_atol

contains

  pure function accs_test_scale(input)

    real(accs_real), intent(in) :: input
    real(accs_real) :: accs_test_scale
    
    accs_test_scale = TESTSCALE * input

  end function accs_test_scale

  pure function accs_test_atol(input)

    real(accs_real), intent(in) :: input
    real(accs_real) :: accs_test_atol

    accs_test_atol = accs_real_eps * input

  end function accs_test_atol

end module accs_test_utils
