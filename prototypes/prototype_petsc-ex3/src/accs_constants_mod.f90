!> @brief Module file accs_constants.mod
!>
!> @details Defines constants for use in ASiMoV-CCS

module accs_constants

  use accs_kinds, only : accs_int

  implicit none

  private
  
  integer(accs_int), public, parameter :: add_mode = 1
  
end module accs_constants
