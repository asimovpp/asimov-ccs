!v Module file globals.mod
!
!  Defines global variables for use in ASiMoV-CCS

module globals


  use types, only: ccs_mesh
  implicit none

  private

  type(ccs_mesh), public :: mesh


end module globals 
