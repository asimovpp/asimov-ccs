!!! -*- mode: F90 -*-
!!! vim: set syntax=fortran:
!!!
!!!        FILE: mod_accs_test_utils.f90
!!! DESCRIPTION: Module providing utility functions/subroutines for the ASiMoV-CCS testing
!!!              framework.

module accs_test_utils

  use constants, only : accs_real, accs_real_eps
  
  implicit none

  private

  public :: accs_test_scale, accs_test_atol

contains

  pure real(accs_real) function accs_test_scale(input)
    !! Function returning a variable scaled to the range of TESTSCALE
    !!
    !! INPUTS:
    !! + input - the input variable, expected |input| <= 1.
    !! RETURNS:
    !! + accs_test_scale - the input scaled to the range of TESTSCALE
    
    real(accs_real), intent(in) :: input
    
    accs_test_scale = TESTSCALE * input

  end function accs_test_scale

  pure real(accs_real) function accs_test_atol(input)
    !! Function returning the absolute tolerance of a variable based on ```accs_real_eps```,
    !! this gives the allowable error magnitude of a real calculation.
    !!
    !! INPUTS:
    !! + input - the value being tested for equality
    !! RETURNS:
    !! + accs_test_atol - the allowable error magnitude.
    
    real(accs_real), intent(in) :: input

    accs_test_atol = accs_real_eps * abs(input)

  end function accs_test_atol

end module accs_test_utils
