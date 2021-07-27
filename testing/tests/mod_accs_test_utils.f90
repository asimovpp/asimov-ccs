!!!        FILE: mod_accs_test_utils.f90
!!! DESCRIPTION: Module providing utility functions/subroutines for the ASiMoV-CCS testing
!!!              framework.

module accs_test_utils

  implicit none

  private

contains

  pure function accs_test_scale(input)

    real(accs_dbl) :: input

    accs_test_scale = TESTSCALE * input

  end function accs_test_scale

end module accs_test_utils
