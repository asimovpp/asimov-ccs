!> Test transient kernel implementation
!
!> The test is based on integrating a known function in time, this allows us to test the transient
!> kernels and their convergence with, and without, startup effects.
!> The function is integrated to T=t0+10dt and the error evaluated based on the known function, the
!> integration is then repeated refining the timestep each time after which the error convergence
!> can be compared against the expected value.
!
!> We test the linear and non-linear cases:
!>     linear: f(t) = sin(at) => f'(t) = a cos(at)
!> non-linear: f(t) = 1 / (t + C) => f'(t) = -f(t)^2; f(t=0) = 1/C /= 0.

program test_transient_kernel

  use testing_lib
  use kinds, only: ccs_real 
  use error_analysis, only: get_order
  use transient_kernel_def, only: transient_second_order_kernel
  use transient_kernels, only: transient_kernel

  implicit none

  integer, parameter :: nref = 4 ! Number of timestep refinements to perform
  integer, parameter :: nstep0 = 10 ! Number of timesteps to perform before refinement
  real(ccs_real), parameter :: alpha = 3.1415_ccs_real ! Arbitrary constant for linear ODE problem
  real(ccs_real), parameter :: C = 1.617_ccs_real ! Arbitrary constant for non-linear ODE problem

  type(transient_second_order_kernel) :: transient ! The transient kernel
  
  real(ccs_real) :: t0, tend ! Start and end of integration interval
  real(ccs_real) :: dt       ! Timestep
  real(ccs_real) :: V        ! Volume
  integer :: nsteps ! Number of timesteps to integrate over
  
  real(ccs_real), dimension(nref + 1) :: err ! The error history
  real(ccs_real), dimension(nref + 1) :: dts ! dt history
  real(ccs_real) :: cur_dt

  real(ccs_real), dimension(:), allocatable :: f0 ! The initial value(s)

  real(ccs_real) :: order
  integer :: i, j

  call init()
  transient = transient_second_order_kernel()
  call transient%set_step(17)

  print *, "order", transient%order
  print *, "width", transient%width
  print *, "implicit", transient%implicit_coeff
  print *, "explicit", transient%explicit_coeffs
  
  ! Allocate space for the old values
  allocate(f0(transient%get_width()))
  
  !! Test a linear ODE
  print *, "Integrating linear ODE f'(t) = a cos(at)"
  dt = 1.0e-3 ! Arbitrary...
  V = dt**(1.0_ccs_real / 3.0_ccs_real) ! Assuming U~=1 would give a CFL=1
  t0 = 0
  tend = t0 + nstep0 * dt

  ! Perform a series of integrations, refining the step each time
  do i = 1, nref + 1
    nsteps = nstep0 * i
    cur_dt = dt/(2**(i-1)) 
    dts(i) = cur_dt

    ! Set initial values
    do j = 1, transient%get_width()
      f0(j) = f1(t0 - (j - 1) * cur_dt)
    end do

    err(i) = integrate(transient, fprime1, fprime_i1, t0, nsteps, cur_dt, f0, f1, 10000, 1e-10_ccs_real)
  end do
  
  ! Check the error convergence
  call get_order(dts, err, order)
  call assert_gt(order, transient%get_order()*0.95_ccs_real, "Convergence order not preserved by transient kernel")

  ! Test a non-linear ODE
  print *, "Integrating non-linear ODE f'(t) = -f(t)^2"
  dt = 1.0e-3 ! Arbitrary...
  V = dt**(1.0_ccs_real / 3.0_ccs_real) ! Assuming U~=1 would give a CFL=1
  t0 = 0
  tend = t0 + nstep0 * dt

  ! Perform a series of integrations, refininng the step each time
  do i = 1, nref + 1
    nsteps = nstep0 * i
    cur_dt = dt/(2**(i-1)) 
    dts(i) = cur_dt

    ! Set initial values
    do j = 1, transient%get_width()
      f0(j) = f2(t0 - (j - 1) * dt)
    end do

    err(i) = integrate(transient, fprime2, fprime_i2, t0, nsteps, cur_dt, f0, f2, 1000, 1e-8_ccs_real)
  end do

  ! Check the error convergence
  call get_order(dts, err, order)
  call assert_gt(order, transient%get_order()*0.95_ccs_real, "Convergence order not preserved by transient kernel")

  call fin()
  
contains

  !! Linear problem
  pure real(ccs_real) function f1(t)
    real(ccs_real), intent(in) :: t

    f1 = sin(alpha * t)
  end function f1
  pure real(ccs_real) function fprime1(f, t)
    real(ccs_real), intent(in) :: f
    real(ccs_real), intent(in) :: t

    associate(foo => f)
    end associate

    fprime1 = alpha * cos(alpha * t)
  end function fprime1
  pure real(ccs_real) function fprime_i1(f, t)
    real(ccs_real), intent(in) :: f
    real(ccs_real), intent(in) :: t

    associate(foo => f, bar => t)
    end associate

    fprime_i1 = 0.0_ccs_real
  end function fprime_i1

  !! Non-linear problem
  pure real(ccs_real) function f2(t)
    real(ccs_real), intent(in) :: t

    f2 = 1.0_ccs_real / (t + C)
  end function f2
  pure real(ccs_real) function fprime2(f, t)
    real(ccs_real), intent(in) :: f
    real(ccs_real), intent(in) :: t

    associate(foo => f, bar => t)
    end associate

    ! Using Picard linearisation -> no explicit forcing term
    fprime2 = 0.0_ccs_real
  end function fprime2
  pure real(ccs_real) function fprime_i2(f, t)
    real(ccs_real), intent(in) :: f
    real(ccs_real), intent(in) :: t

    associate(foo => t)
    end associate

    fprime_i2 = 2 * f
  end function fprime_i2

  !v Integrate the derivative fprime(t) from t0 over nsteps of dt, using old values f0. At the end
  !  of the integration compute the error using fn(t) and return.
  pure real(ccs_real) function integrate(transient, fprime, fprime_i, t0, nsteps, dt, f0, fn, niter, tol)
    class(transient_kernel), intent(in) :: transient
    real(ccs_real), intent(in) :: t0
    integer, intent(in) :: nsteps
    real(ccs_real), intent(in) :: dt
    real(ccs_real), dimension(:), intent(in) :: f0
    integer, intent(in) :: niter      ! How many non-linear iterations?
    real(ccs_real), intent(in) :: tol ! Non-linear tolerance
    interface
      !> Function to evaluate the transient forcing at time t, may be non-linear.
      pure real(ccs_real) function fprime(f, t)
        use kinds, only: ccs_real
        real(ccs_real), intent(in) :: f
        real(ccs_real), intent(in) :: t
      end function fprime
      !> Function to evaluate the implicit component of the forcing at time t, may be non-linear.
      pure real(ccs_real) function fprime_i(f, t)
        use kinds, only: ccs_real
        real(ccs_real), intent(in) :: f
        real(ccs_real), intent(in) :: t
      end function fprime_i
      !> The function being integrated.
      pure real(ccs_real) function fn(t)
        use kinds, only: ccs_real
        real(ccs_real), intent(in) :: t
      end function fn
    end interface

    real(ccs_real) :: f
    real(ccs_real) :: t
    integer :: step
    real(ccs_real), dimension(size(f0)) :: old

    ! Initialise integration
    old = f0
    f = f0(1)

    ! Perform integration
    t = t0
    do step = 1, nsteps
      f = converge_nonlinear(transient, fprime, fprime_i, old, t, dt, niter, tol)
      old = update_old(f, old)

      t = t + dt
    end do
    
    ! Return error
    integrate = abs(f - fn(t))
    
  end function integrate

  pure function update_old(new, old) result(new_old)
    real(ccs_real), intent(in) :: new
    real(ccs_real), dimension(:), intent(in) :: old
    real(ccs_real), dimension(size(old)) :: new_old

    integer :: i

    do i = size(old) - 1, 1, -1
      new_old(i + 1) = old(i)
    end do
    new_old(1) = new
  end function update_old

  pure real(ccs_real) function converge_nonlinear(transient, fprime, fprime_i, old, t, dt, niter, tol)
    class(transient_kernel), intent(in) :: transient ! The transient kernel
    real(ccs_real), dimension(:), intent(in) :: old  ! Function value(s) at previous timestep(s)
    real(ccs_real), intent(in) :: t                  ! Current time
    real(ccs_real), intent(in) :: dt                 ! Timestep
    integer, intent(in) :: niter                     ! Maximum number of non-linear iterations
    real(ccs_real), intent(in) :: tol                ! (Relative) tolerance
    interface
      !> Function to evaluate the transient forcing at time t, may be non-linear.
      pure real(ccs_real) function fprime(f, t)
        use kinds, only: ccs_real
        real(ccs_real), intent(in) :: f
        real(ccs_real), intent(in) :: t
      end function fprime
      !> Function to evaluate the implicit component of the forcing at time t, may be non-linear.
      pure real(ccs_real) function fprime_i(f, t)
        use kinds, only: ccs_real
        real(ccs_real), intent(in) :: f
        real(ccs_real), intent(in) :: t
      end function fprime_i
    end interface

    real(ccs_real) :: prev ! Currently computed function value
    real(ccs_real) :: fnew ! The updated function value
    logical :: converged
    integer :: i           ! Iteration counter
    real(ccs_real), parameter :: rho = 1.0_ccs_real

    real(ccs_real) :: coeff ! The implicit coefficient from transient and forcing contributions
    real(ccs_real) :: rhs   ! The explicit term from transient and forcing contributions

    converged = .false.
    i = 0
    prev = old(1)
    do while((.not. converged) .and. (i < niter))
      ! Get transient coefficient and RHS
      call transient%eval_coeffs(rho, V, dt, coeff)
      call transient%eval_explicit(rho, V, old, dt, rhs)
      
      ! Add forcing coefficient and RHS
      ! Note that these are implicit schemes so should be evaluated at t+dt
      coeff = coeff + fprime_i(prev, t + dt)
      rhs = rhs + fprime(prev, t + dt)

      ! Perform update
      fnew = rhs / coeff

      ! Check convergence
      converged = (abs(fnew - prev) <= abs(tol * prev))
      
      prev = fnew
      i = i + 1
    end do
    
    converge_nonlinear = fnew
    
  end function converge_nonlinear
  
end program test_transient_kernel
