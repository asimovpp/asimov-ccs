!!! -*- mode: F90 -*-
!!! vim: set syntax=fortran:
!!!
  
subroutine calc_adv_coeffs(coeffL, coeffR, rhs, phiL, phiR, advvel)

  use constants
  
  implicit none

  real(accs_real), intent(in) :: phiL, phiR, advvel
  real(accs_real), intent(out) :: coeffL, coeffR, rhs

  coeffL = 1 * phiL
  coeffR = 2 * phiR
  rhs = 3 * advvel
  
end subroutine calc_adv_coeffs
