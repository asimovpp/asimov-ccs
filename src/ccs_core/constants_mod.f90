!> @brief Module file constants.mod
!
!> @details Defines constants for use in ASiMoV-CCS

module constants

  use kinds, only : accs_int

  implicit none

  private
  
  !> @brief Constants to control setting values in objects
  integer(accs_int), public, parameter :: add_mode = 1    !> Add to existing value
  integer(accs_int), public, parameter :: insert_mode = 2 !> Overwrite existing value

end module constants
