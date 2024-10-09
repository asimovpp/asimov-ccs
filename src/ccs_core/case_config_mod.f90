!v Module file case_config.mod
!
!  Provides concrete types and bases of extensible types.

module case_config

  use kinds, only: ccs_int, ccs_real
  use parallel_types, only: parallel_environment

  implicit none

  private

  real(ccs_real), public :: cfl = huge(0.0)

  ! Fluid properties
  real(ccs_real), public :: viscocity = huge(0.0)

  character(len=:), allocatable, save, public :: case_name

  logical, save, public :: write_gradients = .false.

end module case_config
