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

  ! Field types
  integer, public, parameter :: face_centred = 0         !< Indicates face centred variable
  integer, public, parameter :: cell_centred_upwind = 1  !< Indicates cell centred variable (upwind scheme)
  integer, public, parameter :: cell_centred_central = 2 !< Indicates cell centred variable (central scheme)
  integer, public, parameter :: cell_centred_gamma = 3 !< Indicates cell centred variable (gamma scheme)

  ! field names
  integer(ccs_int), public, parameter :: field_u = 0
  integer(ccs_int), public, parameter :: field_v = 1
  integer(ccs_int), public, parameter :: field_w = 2
  integer(ccs_int), public, parameter :: field_p = 3
  integer(ccs_int), public, parameter :: field_p_prime = 4
  integer(ccs_int), public, parameter :: field_mf = 5

  integer(ccs_int), public, parameter :: ccs_string_len = 128

end module constants
