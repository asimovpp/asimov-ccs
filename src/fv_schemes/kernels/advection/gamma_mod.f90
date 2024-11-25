module advection_gamma_mod
  use advection_mod
  use types

  implicit none

  type, extends(advection_kernel) :: gamma_kernel
  contains
    procedure eval_coeffs => advect_gamma_coeffs
    procedure eval_explicit => advect_gamma_eval
    procedure get_width => get_gamma_width
    procedure get_order => get_gamma_order
  end type gamma_kernel

  interface
    module pure function advect_gamma_coeffs(self) result(coeffs)
      class(advection_kernel), intent(in) :: self
      real(ccs_real), allocatable :: coeffs(:)
    end function advect_gamma_coeffs

    module subroutine advect_gamma_eval(self, result)
      class(advection_kernel), intent(in) :: self
      real(ccs_real), intent(out) :: result
    end subroutine advect_gamma_eval

    module pure function get_gamma_width(self) result(width)
      class(advection_kernel), intent(in) :: self
      integer(ccs_int) :: width
    end function get_gamma_width

    module pure function get_gamma_order(self) result(order)
      class(advection_kernel), intent(in) :: self
      integer(ccs_int) :: order
    end function get_gamma_order
  end interface
end module advection_gamma_mod
