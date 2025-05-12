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
    module pure function advection_coeffs(self) result(coeffs)
      class(advection_kernel), intent(in) :: self
      real(ccs_real), allocatable :: coeffs(:)
    end function advection_coeffs

    module subroutine advection_eval(self, result)
      class(advection_kernel), intent(in) :: self
      real(ccs_real), intent(out) :: result
    end subroutine advection_eval

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
    module pure function diffusion_coeffs(self) result(coeffs)
      class(diffusion_kernel), intent(in) :: self
      real(ccs_real), allocatable :: coeffs(:)
    end function diffusion_coeffs

    module subroutine diffusion_eval(self, result)
      class(diffusion_kernel), intent(in) :: self
      real(ccs_real), intent(out) :: result
    end subroutine diffusion_eval

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
