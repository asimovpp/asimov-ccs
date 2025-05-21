submodule(fv_kernels) advection_common
  use kinds, only: ccs_real
  implicit none

contains

  function recast_as_upwind(self) result(concrete_ptr)
    class(advection_kernel), intent(inout) :: self
    type(upwind_advection_kernel), pointer :: concrete_ptr

    select type (self)
    class is (upwind_advection_kernel)
      concrete_ptr => self
    class default
      concrete_ptr => null()
      error stop "Invalid type for advection kernel"
    end select
  end function recast_as_upwind

  !> Base implementation of the advection kernel
  !> Simple prototype
  module procedure advection_eval_coeffs
  coeffs = [1.0_ccs_real, 0.0_ccs_real]  ! Example coefficients
  end procedure advection_eval_coeffs

  module procedure advection_eval_explicit
  result = 0.0_ccs_real  ! Example explicit term
  end procedure advection_eval_explicit

  module procedure advection_width
  width = 1
  end procedure advection_width

  module procedure advection_order
  order = 1
  end procedure advection_order

end submodule advection_mod
