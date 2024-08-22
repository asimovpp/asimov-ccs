module fv_kernels
    use types
    implicit none
  !> Advection kernel
    type, extends(abstract_kernel) :: advection_kernel
    contains
        procedure :: coeffs => advection_coeffs
        procedure :: eval => advection_eval
    end type advection_kernel
end module fv_kernels