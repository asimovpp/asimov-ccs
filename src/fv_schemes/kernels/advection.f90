submodule(fv_kernels) advection_kernel_submodule

   implicit none

contains

   !> Simple prototype
   module procedure advection_coeffs
   coeffs = [1.0, 2.0, 3.0, 4.0]  ! Example coefficients
   end procedure advection_coeffs

   module procedure advection_eval
   result = sum(this%coeffs())
   end procedure advection_eval

end submodule advection_kernel_submodule