!v Module file base_mod
!
!  Defines base variables for use in ASiMoV-CCS

module ccs_base


  use types, only: ccs_mesh
  implicit none

  private

  type(ccs_mesh), public :: mesh


end module ccs_base
