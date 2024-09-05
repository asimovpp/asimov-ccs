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
  end type
  interface transient_first_order_kernel
    module procedure init_first_order
  end interface

  type, extends(transient_kernel) :: transient_second_order_kernel
  end type
  interface transient_second_order_kernel
    module procedure init_second_order
  end interface

  type, extends(transient_kernel) :: transient_theta_kernel
  end type
  interface transient_theta_kernel
    module procedure init_theta
  end interface

  private :: init_first_order
  private :: init_second_order
  private :: init_theta

  contains

    function init_first_order() result(transient)
      type(transient_first_order_kernel) :: transient

      transient%order = 1
      transient%width_trans = [1]
      transient%explicit_coeffs_trans = reshape([1.0_ccs_real], shape=(/1,1/))
      transient%implicit_coeff_trans = [ 1.0_ccs_real ]
   end function


  function init_second_order() result(transient)
    type(transient_second_order_kernel) :: transient

    transient%order = 2
    transient%width_trans = [1, 2]
    transient%explicit_coeffs_trans = reshape([1.0_ccs_real, 0.0_ccs_real, &
                                               2.0_ccs_real, -0.5_ccs_real], shape=(/2, 2/))
    transient%implicit_coeff_trans = [ 1.0_ccs_real, 1.5_ccs_real ]
  end function

  function init_theta() result(transient)
    type(transient_theta_kernel) :: transient
    real(ccs_real), parameter :: theta = 1.5_ccs_real

    transient%order = theta
    transient%width_trans = [1, 2]
    transient%explicit_coeffs_trans = reshape([1.0_ccs_real, 0.0_ccs_real, &
                                               1.0_ccs_real + theta, -0.5_ccs_real*theta], shape=(/2, 2/))
    transient%implicit_coeff_trans = [ 1.0_ccs_real, 1.0_ccs_real + theta]
  end function


end module transient_kernel_def