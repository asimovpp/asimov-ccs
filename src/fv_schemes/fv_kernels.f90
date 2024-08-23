module fv_kernels
  use types
  implicit none
  !> Advection kernel
  type, extends(abstract_kernel) :: advection_kernel
  contains
    procedure :: eval_coeffs => advection_coeffs
    procedure :: eval_explicit => advection_eval
    procedure :: get_width => advection_width
  end type advection_kernel

  interface
    module pure function advection_coeffs(self) result(coeffs)
      class(advection_kernel), intent(in) :: self
      real, allocatable :: coeffs(:)
    end function advection_coeffs

    module subroutine advection_eval(self, result)
      class(advection_kernel), intent(in) :: self
      real, intent(out) :: result
    end subroutine advection_eval

    module pure function advection_width(self) result(width)
      class(advection_kernel), intent(in) :: self
      integer :: width
    end function advection_width
  end interface

end module fv_kernels
