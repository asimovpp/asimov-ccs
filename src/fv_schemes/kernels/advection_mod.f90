module advection_mod
  use types
  use kinds, only: ccs_int, ccs_real, ccs_long
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
    ! Using deferred correction phi_f = (a_P,UD phi_P + a_F,UD phi_F)^{n+1} + {[(a_P - a_P,UD) phi_P + (a_F - a_F,UD) phi_F]}^{n}
    ! So advection_coeffs returns the upwind coefficients for all versions
    ! advection_eval returns the explicit evaluation of the difference between advection scheme and upwind + any correction terms.

    module pure function advection_coeffs(self, mf) result(coeffs)
      class(advection_kernel), intent(in) :: self
      real(ccs_real), dimension(2) :: coeffs
    end function advection_coeffs

    module function advection_eval(self, mf, lf, grads, rvecs, result) ! TODO update return
      class(advection_kernel), intent(in) :: self
      real(ccs_real), intent(out) :: result
    end subroutine advection_eval
    
    module subroutine diff_eval(self, muAdx, grads, rvecs, result)
      class(advection_kernel), intent(in) :: self
      real(ccs_real), intent(out) :: result
    end subroutine diff_eval

    module pure function advection_width(self) result(width)
      class(advection_kernel), intent(in) :: self
      integer(ccs_int) :: width
    end function advection_width

    module pure function advection_order(self) result(order)
      class(advection_kernel), intent(in) :: self
      integer(ccs_int) :: order
    end function
  end interface

end module advection_mod
