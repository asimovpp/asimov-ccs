module transient_kernels

  use types
  use kinds, only: ccs_real, ccs_int
  implicit none

  type, abstract :: transient_kernel
    real(ccs_real) :: order
    integer(ccs_int) :: width
    real(ccs_real), allocatable, dimension(:) :: explicit_coeffs
    real(ccs_real) :: implicit_coeff

    real(ccs_real), allocatable, dimension(:) :: width_trans
    real(ccs_real), allocatable, dimension(:, :) :: explicit_coeffs_trans
    real(ccs_real), allocatable, dimension(:) :: implicit_coeff_trans
  contains
    procedure :: get_order
    procedure :: get_width
    procedure :: set_step
    procedure :: eval_coeffs
    procedure :: eval_explicit
  end type

  contains


  subroutine set_step(self, step, restart)
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


  pure subroutine eval_coeffs(self, rho, V, dt, coeff)
    class(transient_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: rho
    real(ccs_real), intent(in) :: V
    real(ccs_real), intent(in) :: dt
    real(ccs_real), intent(out) :: coeff

    coeff = self%implicit_coeff * rho*V/dt
  end subroutine

  pure subroutine eval_explicit(self, rho, V, old, dt, rhs)
    class(transient_kernel), intent(in) :: self
    real(ccs_real), intent(in) :: rho
    real(ccs_real), intent(in) :: V
    real(ccs_real), dimension(:), intent(in) :: old
    real(ccs_real), intent(in) :: dt
    real(ccs_real), intent(out) :: rhs

    rhs = dot_product(old(1:self%width), self%explicit_coeffs(1:self%width))*rho*V/dt

  end subroutine

  integer(ccs_int) function get_order(self)
    class(transient_kernel) :: self

    get_order = self%order
  end function


  integer(ccs_int) function get_width(self)
    class(transient_kernel) :: self

    get_width = self%width
  end function


end module transient_kernels