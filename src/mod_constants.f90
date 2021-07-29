!!!        FILE: mod_constants.f90
!!! DESCRIPTION: Module defining constants used by ASiMoV-CCS

module constants

  use iso_fortran_env
  
  implicit none

  !! Basic types
  integer, parameter :: accs_real = real64

  !! Useful numbers
  real(accs_real), parameter :: accs_real_eps = 10 * 1.0e-16
  
end module constants
