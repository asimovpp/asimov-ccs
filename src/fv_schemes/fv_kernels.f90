module fv_kernels
#include "ccs_macros.inc"

  use types
  use kinds, only: ccs_real, ccs_int
  use fv_advection_kernels, only: advection_kernel => advection_kernel
  use central_difference_kernel, only: cd_advection_kernel => cd_advection_kernel
  use gamma_kernel, only: gamma_advection_kernel => gamma_advection_kernel
  use linear_upwind_kernel, only: luds_advection_kernel => luds_advection_kernel
  use upwind_kernel, only: upwind_advection_kernel => upwind_advection_kernel
  use utils, only: exit_print

  implicit none

  !> Diffusion kernel
  type, extends(abstract_kernel) :: diffusion_kernel
  contains
    procedure :: eval_coeffs => diffusion_coeffs
    procedure :: eval_explicit => diffusion_eval
    procedure :: get_width => diffusion_width
    procedure :: get_order => diffusion_order
  end type diffusion_kernel

  interface
    module pure function diffusion_coeffs(self, flux_coeff) result(coeffs)
      class(diffusion_kernel), intent(in) :: self
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), dimension(2) :: coeffs
    end function diffusion_coeffs

    module pure function diffusion_eval(self, phi, flux_coeff, lf, rvecs, grads) result(expl)
      class(diffusion_kernel), intent(in) :: self
      real(ccs_real), dimension(2), intent(in) :: phi
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), intent(in) :: lf
      real(ccs_real), dimension(3, 2), intent(in) :: rvecs
      real(ccs_real), dimension(3, 2), intent(in) :: grads
      real(ccs_real) :: expl
    end function diffusion_eval

    module pure function diffusion_width(self) result(width)
      class(diffusion_kernel), intent(in) :: self
      integer(ccs_int) :: width
    end function diffusion_width

    module pure function diffusion_order(self) result(order)
      class(diffusion_kernel), intent(in) :: self
      integer(ccs_int) :: order
    end function diffusion_order
  end interface

  interface
    module function create_advection_kernel(phi) result(kernel)
      class(field), intent(in) :: phi
      class(advection_kernel), allocatable :: kernel
    end function create_advection_kernel
  end interface

contains

  !> Factory returning an advection kernel matching the field discretisation
  module function create_advection_kernel(phi) result(kernel)

    class(field), intent(in) :: phi                        !< Field whose discretisation selects the kernel
    class(advection_kernel), allocatable :: kernel         !< Allocated kernel instance

    select type (phi)
    type is (central_field)
      allocate (cd_advection_kernel :: kernel)
    type is (upwind_field)
      allocate (upwind_advection_kernel :: kernel)
    type is (gamma_field)
      allocate (gamma_advection_kernel :: kernel)
    type is (linear_upwind_field)
      allocate (luds_advection_kernel :: kernel)
    class default
      call error_abort("Invalid velocity field discretisation.")
    end select

  end function create_advection_kernel

end module fv_kernels
