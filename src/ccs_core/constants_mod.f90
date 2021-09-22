!> @brief Module file constants.mod
!!
!! @details Defines constants for use in ASiMoV-CCS

module constants

  use accs_kinds, only : accs_int

  implicit none

  private
  
  !> @brief Constants to control setting values in objects
  !! @var add_mode - add to existing value
  !! @var insert_mode - overwrite existing value
  integer(accs_int), public, parameter :: add_mode = 1, insert_mode = 2 

end module constants
