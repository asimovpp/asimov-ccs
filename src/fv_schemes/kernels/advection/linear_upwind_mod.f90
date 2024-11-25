module advection_luw_mod
  use advection_mod
  use types

  implicit none

  type, extends(advection_kernel) :: luw_kernel
  contains
    procedure eval_coeffs => advect_luw_coeffs
    procedure eval_explicit => advect_luw_eval
    procedure get_width => get_luw_width
    procedure get_order => get_luw_order
  end type luw_kernel
end module advection_luw_mod
