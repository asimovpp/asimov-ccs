module advection_cds_mod
   use advection_mod
   use types

   implicit none

   type, extends(advection_kernel) :: cds_kernel
   contains
      procedure eval_coeffs => advect_cds_coeffs
      procedure eval_explicit => advect_cds_eval
      procedure get_width => get_cds_width
      procedure get_order => get_cds_order
   end type cds_kernel

   interface
      module pure function advect_cds_coeffs(self) result(coeffs)
         class(advection_kernel), intent(in) :: self
         real(ccs_real), allocatable :: coeffs(:)
      end function advect_cds_coeffs

      module subroutine advect_cds_eval(self, result)
         class(advection_kernel), intent(in) :: self
         real(ccs_real), intent(out) :: result
      end subroutine advect_cds_eval

      module pure function get_cds_width(self) result(width)
         class(advection_kernel), intent(in) :: self
         integer(ccs_int) :: width
      end function get_cds_width

      module pure function get_cds_order(self) result(order)
         class(advection_kernel), intent(in) :: self
         integer(ccs_int) :: order
      end function get_cds_order
   end interface

end module advection_cds_mod
