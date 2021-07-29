!!! -*- mode: F90 -*-
!!! vim: set syntax=fortran:
!!!
  
pure subroutine calc_adv_coeffs(coeffL, coeffR, rhs, advvel)

  use constants
  
  implicit none

  real(accs_real), intent(in) :: advvel
  real(accs_real), intent(out) :: coeffL, coeffR, rhs

  coeffL = 0.5 * advvel
  coeffR = coeffL
  rhs = 0
  
end subroutine calc_adv_coeffs
