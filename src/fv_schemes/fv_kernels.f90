module fv_kernels

  use types
  use kinds, only: ccs_real, ccs_int

  implicit none

  !> Advection kernel
  type, extends(abstract_kernel) :: advection_kernel
  contains
    procedure :: eval_coeffs => advection_eval_coeffs
    procedure :: eval_explicit => advection_eval_explicit
    procedure :: get_width => advection_width
    procedure :: get_order => advection_order
  end type advection_kernel

  interface
    module pure function advection_eval_coeffs(self, flux_coeff) result(coeffs)
      class(advection_kernel), intent(in) :: self
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), dimension(2) :: coeffs
    end function advection_eval_coeffs

    module pure function advection_eval_explicit(self, flux_coeff, lf, rvecs, grads, phi_coeffs) result(expl)
      class(advection_kernel), intent(in) :: self
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), intent(in) :: lf
      real(ccs_real), dimension(3, 2), intent(in) :: rvecs
      real(ccs_real), dimension(3, 2), intent(in) :: grads
      real(ccs_real), dimension(2), intent(in), optional :: phi_coeffs
      real(ccs_real) :: expl
    end function advection_eval_explicit

    module pure function advection_width(self) result(width)
      class(advection_kernel), intent(in) :: self
      integer(ccs_int) :: width
    end function advection_width

    module pure function advection_order(self) result(order)
      class(advection_kernel), intent(in) :: self
      integer(ccs_int) :: order
    end function advection_order
  end interface

  !> Upwind Advection Kernel
  !> This kernel is used for the upwind discretisation of the advection term
  !> in the finite volume scheme.
  type, extends(advection_kernel) :: upwind_advection_kernel
  contains
    procedure :: eval_coeffs => advect_upwind_eval_coeffs
    procedure :: eval_explicit => advect_upwind_eval_explicit
    procedure :: get_width => get_upwind_width
    procedure :: get_order => get_upwind_order
  end type upwind_advection_kernel

  interface
    !> Calculates advection coefficient for neighbouring cell using upwind discretisation
    module pure function advect_upwind_eval_coeffs(self, flux_coeff) result(coeffs)
      class(upwind_advection_kernel), intent(in) :: self
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), dimension(2) :: coeffs
    end function advect_upwind_eval_coeffs

    module pure function advect_upwind_eval_explicit(self, flux_coeff, lf, rvecs, grads, phi_coeffs) result(expl)
      class(upwind_advection_kernel), intent(in) :: self
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), intent(in) :: lf
      real(ccs_real), dimension(3, 2), intent(in) :: rvecs
      real(ccs_real), dimension(3, 2), intent(in) :: grads
      real(ccs_real), dimension(2), intent(in), optional :: phi_coeffs
      real(ccs_real) :: expl
    end function advect_upwind_eval_explicit

    module pure function get_upwind_width(self) result(width)
      class(upwind_advection_kernel), intent(in) :: self
      integer(ccs_int) :: width
    end function get_upwind_width

    module pure function get_upwind_order(self) result(order)
      class(upwind_advection_kernel), intent(in) :: self
      integer(ccs_int) :: order
    end function get_upwind_order
  end interface

  !> Central Difference Advection Kernel
  !> This kernel is used for the central difference discretisation of the advection term
  !> in the finite volume scheme.
  type, extends(advection_kernel) :: cd_advection_kernel
    private
    !> Interpolation factor for central difference scheme
    real(ccs_real) :: m_interpol_factor = 0.5_ccs_real  ! sensible default
  contains
    procedure :: eval_coeffs => advect_cd_eval_coeffs
    procedure :: eval_explicit => advect_cd_eval_explicit
    procedure :: get_width => get_cd_width
    procedure :: get_order => get_cd_order
    procedure :: set_interpolation_factor
    procedure :: get_interpolation_factor
  end type cd_advection_kernel

  interface
    !> Calculates advection coefficient for neighbouring cell using central difference discretisation
    module pure function advect_cd_eval_coeffs(self, flux_coeff) result(coeffs)
      class(cd_advection_kernel), intent(in) :: self
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), dimension(2) :: coeffs
    end function advect_cd_eval_coeffs

    module pure function advect_cd_eval_explicit(self, flux_coeff, lf, rvecs, grads, phi_coeffs) result(expl)
      class(cd_advection_kernel), intent(in) :: self
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), intent(in) :: lf
      real(ccs_real), dimension(3, 2), intent(in) :: rvecs
      real(ccs_real), dimension(3, 2), intent(in) :: grads
      real(ccs_real), dimension(2), intent(in), optional :: phi_coeffs

      real(ccs_real) :: expl
    end function advect_cd_eval_explicit

    module pure function get_cd_width(self) result(width)
      class(cd_advection_kernel), intent(in) :: self
      integer(ccs_int) :: width
    end function get_cd_width

    module pure function get_cd_order(self) result(order)
      class(cd_advection_kernel), intent(in) :: self
      integer(ccs_int) :: order
    end function get_cd_order

    module pure subroutine set_interpolation_factor(self, interpol_fact)
      class(cd_advection_kernel), intent(inout) :: self
      real(ccs_real), intent(in) :: interpol_fact
    end subroutine set_interpolation_factor

    module pure function get_interpolation_factor(self) result(interpol_fact)
      class(cd_advection_kernel), intent(in) :: self
      real(ccs_real) :: interpol_fact
    end function get_interpolation_factor

  end interface

  !> Gamma Advection Kernel
  !> This kernel is used for the gamma discretisation of the advection term
  !> in the finite volume scheme.
  type, extends(advection_kernel) :: gamma_advection_kernel
    private
    real(ccs_real) :: beta_m = 0.35_ccs_real
  contains
    procedure :: eval_coeffs => advect_gamma_eval_coeffs
    procedure :: eval_explicit => advect_gamma_eval_explicit
    procedure :: get_width => get_gamma_width
    procedure :: get_order => get_gamma_order
    procedure :: set_beta_m
    procedure :: get_beta_m
  end type gamma_advection_kernel

  interface
    !> Calculates advection coefficient for neighbouring cell using central difference discretisation
    module pure function advect_gamma_eval_coeffs(self, flux_coeff) result(coeffs)
      class(gamma_advection_kernel), intent(in) :: self
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), dimension(2) :: coeffs
    end function advect_gamma_eval_coeffs

    module pure function advect_gamma_eval_explicit(self, flux_coeff, lf, rvecs, grads, phi_coeffs) result(expl)
      class(gamma_advection_kernel), intent(in) :: self
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), intent(in) :: lf
      real(ccs_real), dimension(3, 2), intent(in) :: rvecs
      real(ccs_real), dimension(3, 2), intent(in) :: grads
      real(ccs_real), dimension(2), intent(in), optional :: phi_coeffs

      real(ccs_real) :: expl
    end function advect_gamma_eval_explicit

    module pure function get_gamma_width(self) result(width)
      class(gamma_advection_kernel), intent(in) :: self
      integer(ccs_int) :: width
    end function get_gamma_width

    module pure function get_gamma_order(self) result(order)
      class(gamma_advection_kernel), intent(in) :: self
      integer(ccs_int) :: order
    end function get_gamma_order

    module pure subroutine set_beta_m(self, new_bm)
      class(gamma_advection_kernel), intent(inout) :: self
      real(ccs_real), intent(in) :: new_bm
    end subroutine set_beta_m

    module pure function get_beta_m(self) result(bm)
      class(gamma_advection_kernel), intent(in) :: self
      real(ccs_real) :: bm
    end function get_beta_m

  end interface

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

    module pure function diffusion_eval(self, flux_coeff, lf, rvecs, grads, phi_coeffs) result(expl)
      class(diffusion_kernel), intent(in) :: self
      real(ccs_real), intent(in) :: flux_coeff
      real(ccs_real), intent(in) :: lf
      real(ccs_real), dimension(3, 2), intent(in) :: rvecs
      real(ccs_real), dimension(3, 2), intent(in) :: grads
      real(ccs_real), dimension(2), intent(in), optional :: phi_coeffs
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

end module fv_kernels
