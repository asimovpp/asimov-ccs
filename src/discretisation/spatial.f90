!!! -*- mode: F90 -*-
!!! vim: set syntax=fortran:
!!!
!!!        FILE: spatial.f90
!!! DESCRIPTION: Implements spatial discretisation schemes.
  
pure subroutine calc_adv_coeffs(coeffL, coeffR, advvel)
  !! Given an advecting velocity, compute the advection coefficients for a face using central
  !! differencing.
  !!
  !! INPUTS:
  !! + advvel - the normal velocity at the face
  !! OUTPUTS:
  !! + coeffL - the coefficient for the contribution to the "left" of the face
  !! + coeffR - the coefficient for the contribution to the "right" of the face

  use constants
  
  implicit none

  real(accs_real), intent(in) :: advvel
  real(accs_real), intent(out) :: coeffL, coeffR

  coeffL = 0.5 * advvel
  coeffR = coeffL
  
end subroutine calc_adv_coeffs
