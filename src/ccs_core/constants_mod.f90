!v Module file constants.mod
!
!  Defines constants for use in ASiMoV-CCS

module constants

  use kinds, only: ccs_int

  implicit none

  private

  ! Constants to control setting values in objects
  integer(ccs_int), public, parameter :: add_mode = 1    !< Add to existing value
  integer(ccs_int), public, parameter :: insert_mode = 2 !< Overwrite existing value

  ! String constants
  character(len=4), public, parameter :: geoext = ".geo"
  character(len=18), public, parameter :: adiosconfig = "_adios2_config.xml"
  character(len=12), public, parameter :: ccsconfig = "_config.yaml"

  !> Dimensionality of problems
  integer(ccs_int), public, parameter :: ndim = 3        ! Always 3D

  ! Discretisation schemes
  integer(ccs_int), public, parameter :: uds = 0
  integer(ccs_int), public, parameter :: cds = 1

  ! Location of vector data
  integer(ccs_int), public, parameter :: cell = 1
  integer(ccs_int), public, parameter :: face = 2

  ! Variables
  integer(ccs_int), public, parameter :: VarU = 1  !< u-velocity
  integer(ccs_int), public, parameter :: VarV = 2  !< v-velocity
  integer(ccs_int), public, parameter :: VarW = 3  !< w-velocity
  integer(ccs_int), public, parameter :: VarP = 4  !< pressure

end module constants
