module advection_upwind_mod
  use advection_mod
  use types

  implicit none

  type, extends(advection_kernel) :: upwind_kernel
  contains
    procedure eval_coeffs => advect_upwind_coeffs
    procedure eval_explicit => advect_upwind_eval
    procedure get_width => get_upwind_width
    procedure get_order => get_upwind_order
  end type upwind_kernel
end module advection_upwind_mod
