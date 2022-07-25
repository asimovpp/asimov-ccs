!v Module file bc_constants.mod
!
!  Defines constants for specifying boundary conditions

module bc_constants

  use kinds, only: ccs_int

  implicit none

  integer(ccs_int), parameter :: bc_region_left = -1
  integer(ccs_int), parameter :: bc_region_right = -2
  integer(ccs_int), parameter :: bc_region_bottom = -3
  integer(ccs_int), parameter :: bc_region_top = -4

  integer(ccs_int), parameter :: bc_type_periodic = 1
  integer(ccs_int), parameter :: bc_type_sym = 2
  integer(ccs_int), parameter :: bc_type_dirichlet = 3
  integer(ccs_int), parameter :: bc_type_neumann = 4
  integer(ccs_int), parameter :: bc_type_extrapolate = 5
  integer(ccs_int), parameter :: bc_type_const_grad = 6
  integer(ccs_int), parameter :: bc_type_wall = 7

end module bc_constants
