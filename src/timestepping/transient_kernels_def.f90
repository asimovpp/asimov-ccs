module transient_kernel_def
!> Definition of transient kernels. Each scheme is defined by a theoretical order, width (=number of old values to use) and
!> explicit coefficients associated with the as well as implicit (diagonal) coefficient
!>
!> Because at the start of a computation, fewer 'old' values are available, reverting to a lower order scheme is necessary.
!> This is the reason why each scheme implements their `width`, `explicit_coeffs` and `implicit_coeff` as arrays, 
!> the first element being used for the 1st timestep, etc. until the full scheme can be used.
!>


  use transient_kernels

  type, extends(transient_kernel) :: transient_first_order_kernel
    contains
    procedure :: init => init_first_order
  end type

  type, extends(transient_kernel) :: transient_second_order_kernel
    contains
    procedure :: init => init_second_order
  end type

  type, extends(transient_kernel) :: transient_theta_kernel
    contains
    procedure :: init => init_theta
  end type

  private :: init_first_order
  private :: init_second_order
  private :: init_theta

  contains

  subroutine init_first_order(self)
    class(transient_first_order_kernel) :: self

    self%order = 1
    self%width_trans = [1]
    self%explicit_coeffs_trans = reshape([1.0_ccs_real], shape=(/1,1/))
    self%implicit_coeff_trans = [ 1.0_ccs_real ]
  end subroutine 


  subroutine init_second_order(self)
    class(transient_second_order_kernel) :: self

    self%order = 2
    self%width_trans = [1, 2]
    self%explicit_coeffs_trans = reshape([1.0_ccs_real, 0.0_ccs_real, &
                                          2.0_ccs_real, -0.5_ccs_real], shape=(/2, 2/))
    self%implicit_coeff_trans = [ 1.0_ccs_real, 1.5_ccs_real ]
  end subroutine

  subroutine init_theta(self)
    class(transient_theta_kernel) :: self
    real(ccs_real), parameter :: theta = 1.5_ccs_real

    if (theta == 1.0_ccs_real) then
      self%order = 2
    else
      self%order = 1
    end if
    self%width_trans = [1, 2]
    self%explicit_coeffs_trans = reshape([1.0_ccs_real, 0.0_ccs_real, &
                                          1.0_ccs_real + theta, -0.5_ccs_real*theta], shape=(/2, 2/))
    self%implicit_coeff_trans = [ 1.0_ccs_real, 1.0_ccs_real + theta]
  end subroutine


end module transient_kernel_def
