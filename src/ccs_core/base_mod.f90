!v Module file base_mod
!
!  Defines base variables for use in ASiMoV-CCS

module ccs_base


  use types, only: ccs_mesh
  implicit none

  private

  type(ccs_mesh), public :: mesh

  character(len=128), parameter :: left = "left"
  character(len=128), parameter :: right = "right"
  character(len=128), parameter :: bottom = "bottom"
  character(len=128), parameter :: top = "top"
  character(len=128), parameter :: back = "back"
  character(len=128), parameter :: front = "front"
  character(len=128), dimension(6), parameter, public :: bnd_names_default = &
       [left, right, bottom, top, back, front]
  
end module ccs_base
