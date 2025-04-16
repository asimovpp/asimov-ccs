submodule(fv_kernels) advection_kernel_submodule

   use kinds, only: ccs_real
   implicit none

contains

   !> Simple prototype
   module procedure advection_coeffs
     coeffs = [1.0_ccs_real, 2.0_ccs_real, 3.0_ccs_real, 4.0_ccs_real]  ! Example coefficients
   end procedure advection_coeffs

   module procedure advection_eval
     result = sum(this%coeffs())
   end procedure advection_eval

end submodule advection_kernel_submodule