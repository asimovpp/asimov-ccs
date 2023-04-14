!v Implementation of a first order timestepping scheme with a fixed-size timestep.
submodule(timestepping) timestepping_first_order

  implicit none

  integer(ccs_int), parameter :: num_old_vals = 1 !< the number of old field values the scheme uses

contains

  module subroutine finalise_timestep()

    ! not required for this timestepping scheme

  end subroutine finalise_timestep

  module subroutine reset_timestepping()

    call reset_timestepping_module()

  end subroutine reset_timestepping

  module subroutine apply_timestep(mesh, phi, diag, M, b)

    type(ccs_mesh), intent(in) :: mesh
    class(field), intent(inout) :: phi
    class(ccs_vector), intent(inout) :: diag
    class(ccs_matrix), intent(inout) :: M
    class(ccs_vector), intent(inout) :: b

    if (.not. timestepping_is_active()) then
      return
    end if

    call apply_timestep_first_order(mesh, phi, diag, M, b)

  end subroutine apply_timestep

  module subroutine update_old_values(x)

    use vec, only: get_vector_data, restore_vector_data

    class(field), intent(inout) :: x

    real(ccs_real), dimension(:), pointer :: values_data, old_values_data

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
