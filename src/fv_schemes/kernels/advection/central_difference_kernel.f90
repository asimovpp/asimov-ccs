module central_difference_kernel
  use types
  use kinds, only: ccs_real, ccs_int
  use fv_advection_kernels
  implicit none

  !> Central Difference Advection Kernel
  type, extends(advection_kernel) :: cd_advection_kernel
    private
    !> Interpolation factor for central difference scheme
    real(ccs_real) :: m_interpol_factor = 0.5_ccs_real  ! sensible default
  contains
    procedure :: eval_coeffs => advect_cd_eval_coeffs
    procedure :: eval_explicit => advect_cd_eval_explicit
    procedure :: get_width => get_cd_width
    procedure :: get_order => get_cd_order
    procedure :: set_interpolation_factor
    procedure :: get_interpolation_factor
  end type cd_advection_kernel

  interface
    module pure function advect_cd_eval_coeffs(self, flux_coeff) result(coeffs)
      class(cd_advection_kernel), intent(in) :: self
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), dimension(2) :: coeffs
    end function advect_cd_eval_coeffs

    module pure function advect_cd_eval_explicit(self, phi, flux_coeff, lf, rvecs, grads) result(expl)
      class(cd_advection_kernel), intent(in) :: self
      real(ccs_real), dimension(2), intent(in) :: phi
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), intent(in) :: lf
      real(ccs_real), dimension(3, 2), intent(in) :: rvecs
      real(ccs_real), dimension(3, 2), intent(in) :: grads

      real(ccs_real) :: expl
    end function advect_cd_eval_explicit

    module pure function get_cd_width(self) result(width)
      class(cd_advection_kernel), intent(in) :: self
      integer(ccs_int) :: width
    end function get_cd_width

    module pure function get_cd_order(self) result(order)
      class(cd_advection_kernel), intent(in) :: self
      integer(ccs_int) :: order
    end function get_cd_order

    module pure subroutine set_interpolation_factor(self, interpol_fact)
      class(cd_advection_kernel), intent(inout) :: self
      real(ccs_real), intent(in) :: interpol_fact
    end subroutine set_interpolation_factor

    module pure function get_interpolation_factor(self) result(interpol_fact)
      class(cd_advection_kernel), intent(in) :: self
      real(ccs_real) :: interpol_fact
    end function get_interpolation_factor

  end interface

end module central_difference_kernel
