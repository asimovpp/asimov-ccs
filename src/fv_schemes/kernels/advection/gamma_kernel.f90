module gamma_kernel
  use types
  use kinds, only: ccs_real, ccs_int
  use fv_advection_kernels
  implicit none

  !> Gamma Advection Kernel
  type, extends(advection_kernel) :: gamma_advection_kernel
    private
    real(ccs_real) :: beta_m = 0.35_ccs_real
  contains
    procedure :: eval_coeffs => advect_gamma_eval_coeffs
    procedure :: eval_explicit => advect_gamma_eval_explicit
    procedure :: get_width => get_gamma_width
    procedure :: get_order => get_gamma_order
    procedure :: set_beta_m
    procedure :: get_beta_m
  end type gamma_advection_kernel

  interface
    module pure function advect_gamma_eval_coeffs(self, flux_coeff) result(coeffs)
      class(gamma_advection_kernel), intent(in) :: self
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), dimension(2) :: coeffs
    end function advect_gamma_eval_coeffs

    module pure function advect_gamma_eval_explicit(self, phi, flux_coeff, lf, rvecs, grads) result(expl)
      class(gamma_advection_kernel), intent(in) :: self
      real(ccs_real), dimension(2), intent(in) :: phi
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), intent(in) :: lf
      real(ccs_real), dimension(3, 2), intent(in) :: rvecs
      real(ccs_real), dimension(3, 2), intent(in) :: grads

      real(ccs_real) :: expl
    end function advect_gamma_eval_explicit

    module pure function get_gamma_width(self) result(width)
      class(gamma_advection_kernel), intent(in) :: self
      integer(ccs_int) :: width
    end function get_gamma_width

    module pure function get_gamma_order(self) result(order)
      class(gamma_advection_kernel), intent(in) :: self
      integer(ccs_int) :: order
    end function get_gamma_order

    module pure subroutine set_beta_m(self, new_bm)
      class(gamma_advection_kernel), intent(inout) :: self
      real(ccs_real), intent(in) :: new_bm
    end subroutine set_beta_m

    module pure function get_beta_m(self) result(bm)
      class(gamma_advection_kernel), intent(in) :: self
      real(ccs_real) :: bm
    end function get_beta_m

  end interface

end module gamma_kernel
