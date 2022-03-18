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

  !> @brief String constants
  character(len=4), public, parameter :: geoext = ".geo" 
  character(len=18), public, parameter :: adiosconfig = "_adios2_config.xml" 
  character(len=12), public, parameter :: ccsconfig = "_config.yaml" 

  !> @brief Dimensionality of problems
  integer(accs_int), public, parameter :: ndim = 3        !> Always 3D

  !> @brief Discretisation schemes
  integer(accs_int), public, parameter :: uds = 0
  integer(accs_int), public, parameter :: cds = 1

  !> @brief Location of vector data
  integer(accs_int), public, parameter :: cell = 1
  integer(accs_int), public, parameter :: face = 2
  
end module constants
