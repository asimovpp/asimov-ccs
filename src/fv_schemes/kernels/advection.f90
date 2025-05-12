submodule(fv_kernels) advection_kernel_submodule
  use fv_kernels

  implicit none

contains

  !> Simple prototype
  module procedure advection_coeffs
    coeffs = [1.0, 2.0, 3.0, 4.0]  ! Example coefficients
  end procedure advection_coeffs

  module procedure advection_eval
    result = sum(self%eval_coeffs())
  end procedure advection_eval

  module procedure advection_width
    width = 1
  end procedure advection_width

  module procedure advection_order
    order = 1
  end procedure advection_order

end submodule advection_kernel_submodule
