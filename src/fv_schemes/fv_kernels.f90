module fv_kernels

  use types
  use kinds, only: ccs_real, ccs_int

  implicit none

  !> Advection kernel
  type, extends(abstract_kernel) :: advection_kernel
  contains
    procedure :: eval_coeffs => advection_coeffs
    procedure :: eval_explicit => advection_eval
    procedure :: get_width => advection_width
    procedure :: get_order => advection_order
  end type advection_kernel

  interface
    module pure function advection_coeffs(self, flux_coeff) result(coeffs)
      class(advection_kernel), intent(in) :: self
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), dimension(2) :: coeffs
    end function advection_coeffs

    module pure function advection_eval(self, flux_coeff, lf, rvecs, grads) result(expl)
      class(advection_kernel), intent(in) :: self
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), intent(in) :: lf
      real(ccs_real), dimension(3, 2), intent(in) :: rvecs
      real(ccs_real), dimension(3, 2), intent(in) :: grads
      real(ccs_real):: expl
    end function advection_eval

    module pure function advection_width(self) result(width)
      class(advection_kernel), intent(in) :: self
      integer(ccs_int) :: width
    end function advection_width

    module pure function advection_order(self) result(order)
      class(advection_kernel), intent(in) :: self
      integer(ccs_int) :: order  
    end function advection_order
  end interface

  !> Diffusion kernel
  type, extends(abstract_kernel) :: diffusion_kernel
  contains
    procedure :: eval_coeffs => diffusion_coeffs
    procedure :: eval_explicit => diffusion_eval
    procedure :: get_width => diffusion_width
    procedure :: get_order => diffusion_order
  end type diffusion_kernel

  interface
    module pure function diffusion_coeffs(self, flux_coeff) result(coeffs)
      class(diffusion_kernel), intent(in) :: self
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), dimension(2) :: coeffs
    end function diffusion_coeffs

    module pure function diffusion_eval(self, flux_coeff, lf, rvecs, grads) result(expl)
      class(diffusion_kernel), intent(in) :: self
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), intent(in) :: lf
      real(ccs_real), dimension(3, 2), intent(in) :: rvecs
      real(ccs_real), dimension(3, 2), intent(in) :: grads
      real(ccs_real):: expl
    end function diffusion_eval

    module pure function diffusion_width(self) result(width)
      class(diffusion_kernel), intent(in) :: self
      integer(ccs_int) :: width
    end function diffusion_width

    module pure function diffusion_order(self) result(order)
      class(diffusion_kernel), intent(in) :: self
      integer(ccs_int) :: order  
    end function diffusion_order
  end interface

end module fv_kernels
