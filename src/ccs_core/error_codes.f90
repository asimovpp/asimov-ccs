!v Module file error_codes.mod
!
!  Defines error codes for use in ASiMoV-CCS

module error_codes

  use kinds, only: ccs_int

  implicit none

  private

  ! Error stop codes
  integer, public, parameter :: unknown_scheme = 100
  integer, public, parameter :: unknown_type = 101
  integer, public, parameter :: field_not_found = 102
  integer, public, parameter :: field_index_exceeded = 103 ! Field index exceeds number of flow fields
  integer, public, parameter :: invalid_bc_name = 104 
  integer, public, parameter :: invalid_bc_id = 105 
  integer, public, parameter :: bc_index_not_found = 106 
  integer, public, parameter :: no_access_to_cell = 107 ! Trying to access cell I don't have access to
  integer, public, parameter :: self_not_neighbour = 108 ! Attempt to set self as neighbour
  integer, public, parameter :: invalid_neighbour = 109 ! Neighbour index is not valid
  integer, public, parameter :: unknown_bc_type = 110
  integer, public, parameter :: invalid_component = 111 
  integer, public, parameter :: unknown_mode = 112
  integer, public, parameter :: no_free_entry = 113
  integer, public, parameter :: unknown_bc_name = 114

end module error_codes
