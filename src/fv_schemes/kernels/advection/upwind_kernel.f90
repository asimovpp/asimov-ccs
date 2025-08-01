module upwind_kernel
  use types
  use kinds, only: ccs_real, ccs_int
  use fv_advection_kernels
  implicit none

  !> Upwind Advection Kernel
  type, extends(advection_kernel) :: upwind_advection_kernel
  contains
    procedure :: eval_coeffs => advect_upwind_eval_coeffs
    procedure :: eval_explicit => advect_upwind_eval_explicit
    procedure :: get_width => get_upwind_width
    procedure :: get_order => get_upwind_order
  end type upwind_advection_kernel

  interface
    module pure function advect_upwind_eval_coeffs(self, flux_coeff) result(coeffs)
      class(upwind_advection_kernel), intent(in) :: self
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), dimension(2) :: coeffs
    end function advect_upwind_eval_coeffs

    module pure function advect_upwind_eval_explicit(self, phi, flux_coeff, lf, rvecs, grads) result(expl)
      class(upwind_advection_kernel), intent(in) :: self
      real(ccs_real), dimension(2), intent(in) :: phi
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), intent(in) :: lf
      real(ccs_real), dimension(3, 2), intent(in) :: rvecs
      real(ccs_real), dimension(3, 2), intent(in) :: grads
      real(ccs_real) :: expl
    end function advect_upwind_eval_explicit

    module pure function get_upwind_width(self) result(width)
      class(upwind_advection_kernel), intent(in) :: self
      integer(ccs_int) :: width
    end function get_upwind_width

    module pure function get_upwind_order(self) result(order)
      class(upwind_advection_kernel), intent(in) :: self
      integer(ccs_int) :: order
    end function get_upwind_order
  end interface

end module upwind_kernel
