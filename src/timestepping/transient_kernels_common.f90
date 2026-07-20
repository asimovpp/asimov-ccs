module transient_kernels

  use types
  use kinds, only: ccs_real, ccs_int

  implicit none

  type, abstract :: transient_kernel
    integer(ccs_int) :: order = huge(0_ccs_int)          !< Theoretical order of the scheme
    integer(ccs_int), private :: width = huge(0_ccs_int) !< size of the stencil required (=number of old values)
    real(ccs_real), private :: dt = huge(0.0_ccs_real)    !< time step size
    real(ccs_real), private :: invdt = huge(0.0_ccs_real) !< time step size
    real(ccs_real), private :: implicit_coeff = huge(0.0_ccs_real)        !< lhs/diagonal coefficient
    real(ccs_real), allocatable, dimension(:), private :: explicit_coeffs !< rhs coefficients associated to old values

    ! *_trans variable are set once by the `init` function and define the scheme variables depending on the timestep
    ! This is to revert to a lower width scheme for the first few timesteps when timestep < width
    integer(ccs_int), allocatable, dimension(:) :: width_trans
    real(ccs_real), allocatable, dimension(:, :) :: explicit_coeffs_trans
    real(ccs_real), allocatable, dimension(:) :: implicit_coeff_trans
  contains
    procedure(init), deferred :: init
    procedure :: get_order
    procedure :: get_width
    procedure :: get_dt
    procedure :: set_step
    procedure :: set_dt
    procedure :: eval_coeffs
    procedure :: eval_explicit
    procedure :: cleanup_kernel
  end type

  abstract interface
    subroutine init(self)
      import :: transient_kernel
      class(transient_kernel) :: self
    end subroutine
  end interface


  contains

  subroutine set_step(self, step, restart)
    !! To be run at the begining of every timestep, it sets the right coefficients self%width, explicit_coeffs and implicit_coeff
    class(transient_kernel) :: self
    integer(ccs_int), intent(in) :: step
    logical, optional, intent(in) :: restart
    integer(ccs_int) :: index

    index = min(step, size(self%width_trans))

    if (present(restart)) then
      if (restart) then
        index = size(self%width_trans)
      end if
    end if

    self%width = self%width_trans(index)
    self%explicit_coeffs = self%explicit_coeffs_trans(:, index)
    self%implicit_coeff = self%implicit_coeff_trans(index)

  end subroutine


  pure subroutine eval_coeffs(self, rho, V, coeff)
    !! Computes and returns the implicit coefficient (diagonal coeff)
    class(transient_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: rho
    real(ccs_real), intent(in) :: V
    real(ccs_real), intent(out) :: coeff

    coeff = self%implicit_coeff * rho*V*self%invdt
  end subroutine

  pure subroutine eval_explicit(self, rho, V, old, rhs)
    !! Computes and returns the right hand side (explicit qqt)
    class(transient_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: rho
    real(ccs_real), intent(in) :: V
    real(ccs_real), dimension(:), intent(in) :: old
    real(ccs_real), intent(out) :: rhs

    rhs = dot_product(old(1:self%width), self%explicit_coeffs(1:self%width))*rho*V*self%invdt

  end subroutine

  subroutine set_dt(self, dt)

    !! Setter for  time step size
    class(transient_kernel), intent(inout) :: self
    real(ccs_real), intent(in) :: dt

    if (dt == 0.0_ccs_real) then
      error stop "Time step size is zero" ! Abort to avoid division by zero
    end if
    
    self%dt = dt
    self%invdt = 1.0_ccs_real/dt

  end subroutine

  real(ccs_real) function get_dt(self)
    !! Getter for  time step size
    class(transient_kernel), intent(in) :: self

    get_dt = self%dt
  end function

  integer(ccs_int) function get_order(self)
    !! Getter for the analytical order
    class(transient_kernel), intent(in) :: self

    get_order = self%order
  end function

  integer(ccs_int) function get_width(self)
    !! Getter for the stencil width
    class(transient_kernel), intent(in) :: self

    get_width = self%width
  end function

  subroutine cleanup_kernel(self)
    !! Deallocate allocatable array from kernel object
    class(transient_kernel) :: self

    if (allocated(self%explicit_coeffs)) then
      deallocate(self%explicit_coeffs)
    end if

    if (allocated(self%explicit_coeffs_trans)) then
      deallocate(self%explicit_coeffs_trans)
    end if

    if (allocated(self%implicit_coeff_trans)) then
      deallocate(self%implicit_coeff_trans)
    end if

    if (allocated(self%width_trans)) then
      deallocate(self%width_trans)
    end if

  end subroutine

end module transient_kernels
