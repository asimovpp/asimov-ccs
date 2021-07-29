!!! -*- mode: F90 -*-
!!! vim: set syntax=fortran:
!!!
  
subroutine calc_adv_coeffs(coeffL, coeffR, rhs, phiL, phiR, advvel)

  use constants
  
  implicit none

  real(accs_real), intent(in) :: phiL, phiR, advvel
  real(accs_real), intent(out) :: coeffL, coeffR, rhs

  coeffL = 0.5 * (phiL + phiR) * advvel
  coeffR = coeffL
  rhs = 0
  
end subroutine calc_adv_coeffs
