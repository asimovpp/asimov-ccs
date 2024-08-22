submodule (fv_kernels) advection_kernel_submodule
    use types
    
    implicit none

    contains

    !> Simple prototype
    pure function advection_coeffs(this) result(coeffs)
        class(advection_kernel), intent(in) :: this
        real :: coeffs(:)
        coeffs = [1.0, 2.0, 3.0]  ! Example coefficients
    end function advection_coeffs

    subroutine advection_eval(this, result)
        class(advection_kernel), intent(in) :: this
        real, intent(out) :: result
        result = sum(this%coeffs())
    end subroutine advection_eval

    end submodule advection_kernel_submodule