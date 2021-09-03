!> @brief Module file accs_constants.mod
!>
!> @details Defines constants for use in ASiMoV-CCS

module accs_constants

  use accs_kinds, only : accs_int

  implicit none

  private
  
  !> Constants to control setting values in objects
  !! add_mode: add to existing value
  !! insert_mode: overwrite existing value
  integer(accs_int), public, parameter :: add_mode = 1, insert_mode = 2 

end module accs_constants
