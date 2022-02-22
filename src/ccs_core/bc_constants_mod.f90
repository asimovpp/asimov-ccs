!> @brief Module file bc_constants.mod
!
!> @details Defines constants for specifying boundary conditions

module bc_constants

  use kinds, only: accs_int
  
  implicit none

  integer(accs_int), parameter :: bc_region_left = -1
  integer(accs_int), parameter :: bc_region_right = -2
  integer(accs_int), parameter :: bc_region_top = -3
  integer(accs_int), parameter :: bc_region_bottom = -4

  integer(accs_int), parameter :: bc_type_periodic = 1
  integer(accs_int), parameter :: bc_type_sym = 2
  integer(accs_int), parameter :: bc_type_dirichlet = 3
  integer(accs_int), parameter :: bc_type_const_grad = 4

end module bc_constants
