!v Implementation of a first order timestepping scheme with a fixed-size timestep.
submodule(timestepping) timestepping_first_order

  implicit none

  integer(ccs_int), parameter :: num_old_vals = 1 !< the number of old field values the scheme uses
  real(ccs_real), parameter :: theoretical_order = 1.0_ccs_real !< Theoretical order of accuracy of the scheme

contains

  module subroutine finalise_timestep()

    call increment_time_step()

  end subroutine finalise_timestep

  module subroutine reset_timestepping()

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

    if (.not. timestepping_is_active()) then
      return
    end if

    call apply_timestep_first_order(phi, diag, M, b)

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
