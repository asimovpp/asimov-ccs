!v Implementation of the theta timestepping scheme with a fixed-size timestep.
submodule(timestepping) timestepping_theta

  use transient_kernel_def, only: transient_theta_kernel
  implicit none

  real(ccs_real), parameter :: theta = 0.5 !< timestepping scheme mixing factor
  integer(ccs_int), parameter :: num_old_vals = 2 !< the number of old field values the scheme uses
  real(ccs_real), parameter :: theoretical_order = theta + 1.0_ccs_real !< Theoretical order of accuracy of the scheme
  logical, save :: first_update = .true.

contains

  module subroutine finalise_timestep()

    first_update = .false.
    call increment_time_step()

  end subroutine finalise_timestep

  module subroutine reset_timestepping()

    first_update = .true.
    call reset_timestepping_module()

  end subroutine reset_timestepping

  pure module subroutine get_theoretical_order(order)
    real(ccs_real), intent(out) :: order

    order = theoretical_order

  end subroutine

  module subroutine apply_timestep(phi, diag, M, b)

    class(field), intent(inout) :: phi
    class(ccs_vector), intent(inout) :: diag
    class(ccs_matrix), intent(inout) :: M
    class(ccs_vector), intent(inout) :: b
    class(transient_kernel), allocatable :: transient ! The transient kernel

    if (.not. timestepping_is_active()) then
      return
    end if

    allocate(transient_theta_kernel :: transient)
    call apply_timestep_kernel(transient, phi, diag, M, b)

  end subroutine apply_timestep

  module subroutine update_old_values(x)

    class(field), intent(inout) :: x

    if (.not. timestepping_is_active()) then
      return
    end if

    call update_old_values_generic(num_old_vals, x)

  end subroutine

  module subroutine initialise_old_values(vec_properties, x)

    type(vector_spec), intent(in) :: vec_properties
    class(field), intent(inout) :: x

    call initialise_old_values_generic(vec_properties, num_old_vals, x)

  end subroutine

end submodule

