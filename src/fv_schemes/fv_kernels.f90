module fv_kernels
  use types
  use kinds, only: ccs_real
  implicit none
  !> Advection kernel
  type, extends(abstract_kernel) :: advection_kernel
  contains
     procedure :: coeffs => advection_coeffs
     procedure :: eval => advection_eval
  end type advection_kernel

  interface
     module pure function advection_coeffs(this) result(coeffs)
        class(advection_kernel), intent(in) :: this
        real(ccs_real), allocatable :: coeffs(:)
     end function advection_coeffs

     module subroutine advection_eval(this, result)
        class(advection_kernel), intent(in) :: this
        real(ccs_real), intent(out) :: result
     end subroutine advection_eval
  end interface

end module fv_kernels