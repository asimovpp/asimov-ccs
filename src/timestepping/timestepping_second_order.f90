!v Implementation of a second order timestepping scheme with a fixed-size timestep.
submodule(timestepping) timestepping_second_order

  use transient_kernel_def, only: transient_second_order_kernel
  implicit none

  integer(ccs_int), parameter :: num_old_vals = 2 !< the number of old field values the scheme uses
  logical, save :: first_update = .true.
  real(ccs_real), parameter :: theoretical_order = 2.0_ccs_real !< Theoretical order of accuracy of the scheme

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

    ! Allocating the kernel every timestep is temporary, it should be done once when the final equation 
    ! system isget_vector_data, restore_vector_data implemented
    allocate(transient_second_order_kernel :: transient)
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
