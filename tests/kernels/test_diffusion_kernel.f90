!> Test diffusion kernel implementation
!
!> The test is based on interpolating a known function over an interval, this allows us to test the diffusion
!> kernels.

program test_diffusion_kernel

    use kinds, only: ccs_real
  
    implicit none
  
    integer, parameter :: nref = 4 ! Number of interval refinements to perform
    integer, parameter :: nstep0 = 10 ! Number of timesteps to perform before refinement
    real(ccs_real), parameter :: alpha = 3.1415_ccs_real ! Arbitrary constant for linear ODE problem
    real(ccs_real), parameter :: C = 1.617_ccs_real ! Arbitrary constant for non-linear ODE problem
  
    type(diffusion_kernel) :: diffusion ! The diffusion kernel
    
    real(ccs_real) :: x0, xend ! Start and end of interval
    real(ccs_real) :: dx       ! mesh interval
    integer :: nsteps ! Number of mesh intervals to integrate over
    
    real(ccs_real) :: ctol, r
    real(ccs_real), dimension(nref + 1) :: err ! The error history
  
    real(ccs_real), dimension(:), allocatable :: f0 ! The initial value(s)
  
    integer :: i
  
    ! The convergence ratio between refinement levels should be (1/2)^p where p is the order of the
    ! scheme. To allow some reasonable tolerance, 1.1/(2^p) has been used.
    ctol = 1.1 / (2**diffusion%get_order())
  
    ! Allocate space for the old values
    allocate(f0(diffusion%get_width()))
    
  contains
  
    !! Linear problem
    pure real(ccs_real) function f1(x)
      real(ccs_real), intent(in) :: x
  
      f1 = alpha * x + C
    end function f1
    
  end program test_diffusion_kernel
  
  