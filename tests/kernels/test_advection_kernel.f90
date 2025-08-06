program test_advection_kernels
  use testing_lib
  use kinds,          only : ccs_real, ccs_int
  use utils,          only : str
  use error_analysis, only : compute_order
  use fv_kernels

  implicit none

  !---------------------------- configuration --------------------------------
  integer,        parameter :: nLevels = 6           ! refinements (Δx = L/2^ℓ)
  real(ccs_real), parameter :: L       = 1.0_ccs_real
  real(ccs_real), parameter :: u0      = 2.0_ccs_real  ! constant face velocity
  real(ccs_real), parameter :: phi0    = 3.14_ccs_real ! for const‑φ test
  real(ccs_real), parameter :: interpF = 0.5_ccs_real  ! generic λ_f

  ! kernel instances
  type(gamma_advection_kernel)   :: gamma
  type(upwind_advection_kernel)  :: upwind
  type(luds_advection_kernel)    :: luds
  type(cd_advection_kernel)      :: cd

  ! error containers
  real(ccs_real) :: h(nLevels),  &
                    err_gamma(nLevels), err_upwind(nLevels), &
                    err_luds(nLevels),  err_cd(nLevels)

  ! local scratch
  integer(ccs_int) :: lev
  real(ccs_real)   :: dx
  real(ccs_real),  dimension(3) :: xP, xF, xFace
  real(ccs_real)   :: phiP, phiF, phiFace, gradFace
  real(ccs_real),  dimension(3) :: gradP, gradF
  real(ccs_real),  dimension(3,2) :: rvecs, grads
  real(ccs_real),  dimension(2)   :: phiCell
  character(len=*), parameter :: hdr = 'Advection‑kernel test ‖ '

  ! ========================== ACCURACY CHECK =============================
  do lev = 1, nLevels
     dx     = L / real(2**lev, ccs_real)
     h(lev) = dx

     ! cell centres and face --------------------------------------------------
     xP = [0.0_ccs_real, 0.0_ccs_real, 0.0_ccs_real]
     xF = [dx,          0.0_ccs_real, 0.0_ccs_real]
     xFace = [dx/2.0_ccs_real, 0.0_ccs_real, 0.0_ccs_real]

     ! analytic scalar and gradient ------------------------------------------
     call sample_field(xP(1), phiP, gradP(1))
     call sample_field(xF(1), phiF, gradF(1))
     call sample_field(xFace(1), phiFace, gradFace)
     gradP(2:3) = 0.0_ccs_real;  gradF(2:3) = 0.0_ccs_real

     ! geometry from face to centres (consistent with kernel impl.)
     rvecs(:,1) = xP - xFace         ! face → owner P
     rvecs(:,2) = xF - xFace         ! face → neighbour F

     grads(:,1) = gradP           ! ∇φ_P
     grads(:,2) = gradF           ! ∇φ_F

     phiCell = [phiP, phiF]       ! [P , F]

     ! ---------------- kernels ----------------------------------------------
     call kernel_error(gamma , phiCell, phiFace, err_gamma(lev))
     call kernel_error(upwind, phiCell, phiFace, err_upwind(lev))
     call kernel_error(luds  , phiCell, phiFace, err_luds (lev))
     call kernel_error(cd    , phiCell, phiFace, err_cd   (lev))
  end do

  call check_order('Gamma'     , h(2:), err_gamma(2:), gamma %get_order())
  call check_order('Upwind', h(2:), err_upwind(2:), upwind%get_order())
  call check_order('LUDS'  , h(2:), err_luds (2:), luds  %get_order())
  call check_order('CD'    , h(2:), err_cd   (2:), cd    %get_order())

  ! =========================== CONSERVATION ===============================
  call check_constant_phi(gamma , 'Gamma')
  call check_constant_phi(upwind, 'Upwind')
  call check_constant_phi(luds  , 'LUDS')
  call check_constant_phi(cd    , 'CD')

  ! ==================== FLOW‑DIRECTION SYMMETRY ===========================
  call check_reverse_flow(gamma , 'Gamma')
  call check_reverse_flow(upwind, 'Upwind')
  call check_reverse_flow(luds  , 'LUDS')
  call check_reverse_flow(cd    , 'CD')


contains
  pure subroutine sample_field(x, phi, dphidx)
    real(ccs_real), intent(in)  :: x
    real(ccs_real), intent(out) :: phi, dphidx
    real(ccs_real), parameter   :: twopi = 2.0_ccs_real * acos(-1.0_ccs_real)
    phi    = sin(twopi*x + acos(-1.0_ccs_real)/7.0_ccs_real)
    dphidx = twopi * cos(twopi*x + acos(-1.0_ccs_real)/7.0_ccs_real)
  end subroutine sample_field

  subroutine kernel_error(k, phiCell, phiFace, err)
    class(advection_kernel), intent(in) :: k
    real(ccs_real), dimension(2), intent(in) :: phiCell   ! [P , F]
    real(ccs_real),               intent(in) :: phiFace   ! analytic φ at face
    real(ccs_real),               intent(out):: err
    real(ccs_real)                            :: flux_disc, flux_exact
    real(ccs_real), dimension(2)              :: coeffs

    coeffs     = k%eval_coeffs(u0)
    flux_disc  = coeffs(2)*phiCell(2) + coeffs(1)*phiCell(1) + &
                 k%eval_explicit(phiCell, u0, interpF, rvecs, grads)
    flux_exact = u0 * phiFace
    err        = abs(flux_exact - flux_disc)
  end subroutine kernel_error

  subroutine check_order(name, h, errors, theo)
    character(*), intent(in) :: name
    real(ccs_real), intent(in) :: h(:), errors(:)
    integer(ccs_int), intent(in) :: theo
    real(ccs_real) :: p, thresh
    logical :: ok

    call compute_order(h, errors, p)
    thresh = 0.85_ccs_real * real(theo,ccs_real)
    call assert_ge(p, thresh, hdr//trim(name)//': order '//str(p)//' < '//str(thresh), outval=ok)
    if (.not. ok) write(*,*) '  ORDER  failure - ', trim(name), ': p =', p
  end subroutine check_order

  subroutine check_constant_phi(k, name)
    class(advection_kernel), intent(in) :: k
    character(*),           intent(in) :: name
    real(ccs_real), dimension(2) :: coeffs
    real(ccs_real) :: rhs, flux
    logical :: ok
    real(ccs_real), dimension(2) :: phi_const
    real(ccs_real), dimension(3,2) :: grads_zero


    coeffs   = k%eval_coeffs(u0)
    phi_const = [phi0, phi0]
    grads_zero = 0.0_ccs_real
    rhs      = k%eval_explicit(phi_const, u0, interpF, rvecs, grads_zero)
    flux     = coeffs(1)*phi0 + coeffs(2)*phi0 + rhs

    call assert_eq(flux, u0*phi0, hdr//trim(name)//': constant-phi flux', outval=ok)
    if (.not. ok) write(*,*) '  CONSERVATION failure - ', trim(name)
  end subroutine check_constant_phi

  subroutine check_reverse_flow(k, name)
    class(advection_kernel), intent(in) :: k
    character(*),           intent(in) :: name
    real(ccs_real), dimension(2) :: coeffs
    real(ccs_real) :: rhs, flux_pos, flux_neg
    logical :: ok
    real(ccs_real), dimension(2) :: phi_const
    real(ccs_real), dimension(3,2) :: grads_zero

    phi_const  = [phi0, phi0]
    grads_zero = 0.0_ccs_real

    ! positive velocity
    coeffs   = k%eval_coeffs(u0)
    rhs      = k%eval_explicit(phi_const, u0, interpF, rvecs, grads_zero)
    flux_pos = coeffs(1)*phi0 + coeffs(2)*phi0 + rhs

    ! negative velocity
    coeffs   = k%eval_coeffs(-u0)
    rhs      = k%eval_explicit(phi_const, -u0, interpF, rvecs, grads_zero)
    flux_neg = coeffs(1)*phi0 + coeffs(2)*phi0 + rhs

    call assert_eq(flux_neg, -flux_pos, hdr//trim(name)//': antisymmetry', outval=ok)
    if (.not. ok) write(*,*) '  ANTISYMMETRY failure - ', trim(name)
  end subroutine check_reverse_flow

end program test_advection_kernels
