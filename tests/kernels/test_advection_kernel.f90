program test_advection_kernels
  use testing_lib
  use kinds,          only : ccs_real, ccs_int
  use utils,          only : str
  use error_analysis, only : compute_order
  use fv_kernels

  implicit none

  ! -----------------------------------------------------------------------
  !  Numerical‑experiment parameters
  ! -----------------------------------------------------------------------
  integer,            parameter :: nLevels = 6      ! number of mesh refinements
  real(ccs_real),     parameter :: L       = 1.0_ccs_real
  real(ccs_real),     parameter :: u0      = 2.0_ccs_real   ! constant face velocity
  real(ccs_real),     parameter :: phi0    = 3.14_ccs_real  ! for conservation tests
  real(ccs_real),     parameter :: interpF = 0.5_ccs_real   ! not used by all schemes

  !  Kernel instances
  type(gamma_advection_kernel)   :: gamma
  type(upwind_advection_kernel)  :: upwind
  type(luds_advection_kernel)    :: luds
  type(cd_advection_kernel)      :: cd

  !  Storage for grid size and errors per refinement level
  real(ccs_real), dimension(nLevels) :: h,   &
                                        err_gamma, err_upwind, &
                                        err_luds,  err_cd

  !  Local scratch variables
  integer(ccs_int) :: lev
  real(ccs_real)   :: dx, phiP, phiN, phiF
  real(ccs_real),   dimension(2)   :: coeffs
  real(ccs_real),   dimension(3,2) :: rvecs, grads
  real(ccs_real),   dimension(3)   :: xP, xN, xf
  real(ccs_real)   :: order
  character(len=*), parameter :: msgHdr = 'Advection‑kernel test ‖ '

  !  r‑vectors / gradients are unused in these 1‑D tests
  rvecs = 0.0_ccs_real;  grads = 0.0_ccs_real

  call init()   ! testing_lib initialisation

  ! *********************************************************************
  !  1.  Convergence‑rate test (uniform grid refinement Δx = L / 2^ℓ)
  ! *********************************************************************
  do lev = 1, nLevels
     dx      = L / (2.0_ccs_real**lev)
     h(lev)  = dx

     xP      = [0.0_ccs_real, 0.0_ccs_real, 0.0_ccs_real]
     xN      = xP + [dx, 0.0_ccs_real, 0.0_ccs_real]
     xf      = xP + [dx/2.0_ccs_real, 0.0_ccs_real, 0.0_ccs_real]

     !  Analytic field   φ(x) = sin(2πx + π/7)
     phiP    = sin(2.0_ccs_real*acos(-1.0_ccs_real)*xP(1) + acos(-1.0_ccs_real)/7.0_ccs_real)
     phiN    = sin(2.0_ccs_real*acos(-1.0_ccs_real)*xN(1) + acos(-1.0_ccs_real)/7.0_ccs_real)
     phiF    = sin(2.0_ccs_real*acos(-1.0_ccs_real)*xf(1) + acos(-1.0_ccs_real)/7.0_ccs_real)

     !  γ‑scheme
     coeffs  = gamma %eval_coeffs(u0)
     call record_error(coeffs, gamma%eval_explicit(u0, interpF, rvecs, grads), &
                       u0, phiP, phiN, phiF, err_gamma(lev))

     !  Upwind
     coeffs  = upwind%eval_coeffs(u0)
     call record_error(coeffs, upwind%eval_explicit(u0, interpF, rvecs, grads), &
                       u0, phiP, phiN, phiF, err_upwind(lev))

     !  LUDS
     coeffs  = luds %eval_coeffs(u0)
     call record_error(coeffs, luds%eval_explicit(u0, interpF, rvecs, grads), &
                       u0, phiP, phiN, phiF, err_luds(lev))

     !  Central‑difference
     coeffs  = cd    %eval_coeffs(u0)
     call record_error(coeffs, cd   %eval_explicit(u0, interpF, rvecs, grads), &
                       u0, phiP, phiN, phiF, err_cd(lev))
  end do

  ! ------ evaluate orders (skip the coarsest level) -------------------
  call check_order('γ'     , h(2:), err_gamma(2:), gamma %get_order())
  call check_order('Upwind', h(2:), err_upwind(2:), upwind%get_order())
  call check_order('LUDS'  , h(2:), err_luds (2:), luds  %get_order())
  call check_order('CD'    , h(2:), err_cd   (2:), cd    %get_order())

  ! *********************************************************************
  !  2.  Conservation (φ ≡ const.)
  ! *********************************************************************
  call check_constant_phi(gamma , 'γ')
  call check_constant_phi(upwind, 'Upwind')
  call check_constant_phi(luds  , 'LUDS')
  call check_constant_phi(cd    , 'CD')

  ! *********************************************************************
  !  3.  Flow‑direction symmetry
  ! *********************************************************************
  call check_reverse_flow(gamma , 'γ')
  call check_reverse_flow(upwind, 'Upwind')
  call check_reverse_flow(luds  , 'LUDS')
  call check_reverse_flow(cd    , 'CD')

contains

  !> -------------------------------------------------------------------
  !>  Discrete vs analytic flux error
  !> -------------------------------------------------------------------
  subroutine record_error(coeffs, rhs, uf, phiP, phiN, phiF, err)
     real(ccs_real), intent(in)  :: coeffs(2), rhs, uf, phiP, phiN, phiF
     real(ccs_real), intent(out) :: err
     real(ccs_real)              :: flux_exact

     flux_exact = uf * phiF
     err        = abs(flux_exact - (coeffs(1)*phiP + coeffs(2)*phiN + rhs))
  end subroutine record_error

  !> -------------------------------------------------------------------
  !>  Observed convergence order must be ≥ 70 % of theoretical.
  !> -------------------------------------------------------------------
  subroutine check_order(name, h, errors, theo)
     character(*),    intent(in) :: name
     real(ccs_real),  intent(in) :: h(:), errors(:)
     integer(ccs_int),intent(in) :: theo
     real(ccs_real)              :: p, lower

     call compute_order(h, errors, p)
     lower = 0.70_ccs_real * real(theo, ccs_real)
     call assert_ge(p, lower, msgHdr // trim(name) //                  &
          ': observed order p=' // str(p) // ' lower than expected')
  end subroutine check_order

  !> -------------------------------------------------------------------
  !>  Constant‑φ conservation check
  !> -------------------------------------------------------------------
  subroutine check_constant_phi(k, name)
     class(advection_kernel), intent(in) :: k
     character(*),           intent(in) :: name
     real(ccs_real), dimension(2) :: coeffs
     real(ccs_real) :: rhs, flux

     coeffs = k%eval_coeffs(u0)
     rhs    = k%eval_explicit(u0, interpF, rvecs, grads)
     flux   = coeffs(1)*phi0 + coeffs(2)*phi0 + rhs

     call assert_eq(flux, u0*phi0, &
          msgHdr // trim(name) // ': φ ≡ const – non‑conservative')
  end subroutine check_constant_phi

  !> -------------------------------------------------------------------
  !>  Flux antisymmetry with respect to flow direction
  !> -------------------------------------------------------------------
  subroutine check_reverse_flow(k, name)
     class(advection_kernel), intent(in) :: k
     character(*),           intent(in) :: name
     real(ccs_real), dimension(2) :: coeffs
     real(ccs_real) :: rhs, flux_pos, flux_neg

     coeffs   = k%eval_coeffs(u0)
     rhs      = k%eval_explicit(u0, interpF, rvecs, grads)
     flux_pos = coeffs(1)*phi0 + coeffs(2)*phi0 + rhs

     coeffs   = k%eval_coeffs(-u0)
     rhs      = k%eval_explicit(-u0, interpF, rvecs, grads)
     flux_neg = coeffs(1)*phi0 + coeffs(2)*phi0 + rhs

     call assert_eq(flux_neg, -flux_pos, &
          msgHdr // trim(name) // ': F(-u_f) ≠ -F(u_f)')
  end subroutine check_reverse_flow

end program test_advection_kernels
