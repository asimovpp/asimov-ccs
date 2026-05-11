!v Module file bc_constants.mod
!
!  Defines constants for specifying boundary conditions

module bc_constants

  use kinds, only: ccs_int

  implicit none

  integer(ccs_int), parameter :: bc_type_dirichlet = 1
  integer(ccs_int), parameter :: bc_type_neumann = 2
  integer(ccs_int), parameter :: bc_type_extrapolate = 3
  integer(ccs_int), parameter :: bc_type_sym = 4
  integer(ccs_int), parameter :: bc_type_inflow = 5
  integer(ccs_int), parameter :: bc_type_outflow = 6
  integer(ccs_int), parameter :: bc_type_wall = 7
  integer(ccs_int), parameter :: bc_type_slip_wall = 8
  integer(ccs_int), parameter :: bc_type_periodic = 9
  integer(ccs_int), parameter :: bc_type_profile = 10

end module bc_constants
