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




end module error_codes
