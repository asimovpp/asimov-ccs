!v Implementation of a second order timestepping scheme with a fixed-size timestep.
submodule(timestepping) timestepping_second_order

  implicit none

  integer(ccs_int), parameter :: num_old_vals = 2 !< the number of old field values the scheme uses
  logical, save :: first_update = .true.

contains

  module subroutine finalise_timestep()

    first_update = .false.

  end subroutine finalise_timestep

  module subroutine apply_timestep(mesh, phi, diag, M, b)

    type(ccs_mesh), intent(in) :: mesh
    class(field), intent(inout) :: phi
    class(ccs_vector), intent(inout) :: diag
    class(ccs_matrix), intent(inout) :: M
    class(ccs_vector), intent(inout) :: b

    if (.not. timestepping_is_active()) then
      return
    end if

    ! Do a first order update the first time because there no two past time steps yet.
    if (first_update) then
      call apply_timestep_first_order(mesh, phi, diag, M, b)
    else
      call apply_timestep_second_order(mesh, phi, diag, M, b)
    end if

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
