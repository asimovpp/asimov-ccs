submodule(timestepping) timestepping_common
#include "ccs_macros.inc"

  use meshing, only: get_local_num_cells, create_cell_locator, get_volume
  use types, only: cell_locator
  use utils, only: exit_print

  implicit none

  logical :: timestepping_active = .false. !< flag to signify whether timestepping should occur
  logical :: timestep_is_set = .false. !< flag to signify whether dt has already been set
  real(ccs_real) :: dt !< timestep size
  integer(ccs_int) :: current_step = 0

contains

  module subroutine activate_timestepping()
    timestepping_active = .true.
  end subroutine

  pure module function timestepping_is_active() result(active)
    logical :: active
    active = timestepping_active
  end function

  module subroutine reset_timestepping_module()

    timestepping_active = .false.
    timestep_is_set = .false.
    current_step = 0

  end subroutine

  module subroutine set_timestep(timestep)

    real(ccs_real), intent(in) :: timestep

    if (.not. timestep_is_set) then
      dt = timestep
      timestep_is_set = .true.
    else
      call error_abort("Attempted to change timestep after it had already been set.")
    end if

  end subroutine set_timestep

  module function get_timestep() result(timestep)

    real(ccs_real) :: timestep

    if (timestep_is_set) then
      timestep = dt
    else
      call error_abort("Attempted to retrieve timestep before it has been set.")
      timestep = -1
    end if

  end function

  pure module subroutine get_current_step(step)

    integer(ccs_int), intent(out) :: step

    if (timestepping_active) then
      step = current_step
    else
      step = -1
    end if

  end subroutine

  pure module subroutine get_current_time(time)

    real(ccs_real), intent(out) :: time

    if (timestepping_active .and. timestep_is_set) then
      time = current_step * dt
    else
      time = -1.0_ccs_real
    end if

  end subroutine

  module subroutine increment_time_step()

    current_step = current_step + 1

  end subroutine

  module subroutine initialise_old_values_generic(vec_properties, num_old_vals, x)

    use types, only: vector_spec
    use vec, only: create_vector
    use utils, only: update

    type(vector_spec), intent(in) :: vec_properties
    integer(ccs_int), intent(in) :: num_old_vals
    class(field), intent(inout) :: x

    integer(ccs_int) :: i

    if (.not. allocated(x%old_values)) then
      allocate (x%old_values(num_old_vals))
    end if

    do i = 1, num_old_vals
      call create_vector(vec_properties, x%old_values(i)%vec)
      call update(x%old_values(i)%vec)
    end do

  end subroutine

  module subroutine update_old_values_generic(num_old_vals, x)

    use vec, only: get_vector_data, restore_vector_data

    integer(ccs_int), intent(in) :: num_old_vals
    class(field), intent(inout) :: x

    real(ccs_real), dimension(:), pointer :: values_data, old_values_data
    integer(ccs_int) :: i

    do i = num_old_vals, 2, -1
      call get_vector_data(x%old_values(i)%vec, old_values_data)
      call get_vector_data(x%old_values(i - 1)%vec, values_data)
      old_values_data = values_data
      call restore_vector_data(x%old_values(i)%vec, old_values_data)
      call restore_vector_data(x%old_values(i - 1)%vec, values_data)
    end do

    call get_vector_data(x%old_values(1)%vec, old_values_data)
    call get_vector_data(x%values, values_data)
    old_values_data = values_data
    call restore_vector_data(x%old_values(1)%vec, old_values_data)
    call restore_vector_data(x%values, values_data)

  end subroutine

  module subroutine apply_timestep_first_order(phi, diag, M, b)

    use kinds, only: ccs_int
    use mat, only: set_matrix_diagonal, get_matrix_diagonal
    use vec, only: get_vector_data, restore_vector_data
    use utils, only: update, finalise

    class(field), intent(inout) :: phi
    class(ccs_vector), intent(inout) :: diag
    class(ccs_matrix), intent(inout) :: M
    class(ccs_vector), intent(inout) :: b

    real(ccs_real), dimension(:), pointer :: diag_data
    real(ccs_real), dimension(:), pointer :: b_data
    real(ccs_real), dimension(:), pointer :: phi_data
    integer(ccs_int) :: i
    integer(ccs_int) :: local_num_cells

    real(ccs_real) :: V_p
    type(cell_locator) :: loc_p

    call finalise(M)
    call get_matrix_diagonal(M, diag)

    call get_vector_data(phi%old_values(1)%vec, phi_data)
    call get_vector_data(diag, diag_data)
    call update(b)
    call get_vector_data(b, b_data)

    call get_local_num_cells(local_num_cells)
    do i = 1, local_num_cells
      call create_cell_locator(i, loc_p)
      call get_volume(loc_p, V_p)

      ! A = A + V/dt
      diag_data(i) = diag_data(i) + V_p / get_timestep()

      ! b = b + V/dt * phi_old
      b_data(i) = b_data(i) + V_p / get_timestep() * phi_data(i)
    end do
    call restore_vector_data(phi%old_values(1)%vec, phi_data)
    call restore_vector_data(diag, diag_data)
    call restore_vector_data(b, b_data)
    call set_matrix_diagonal(diag, M)

  end subroutine apply_timestep_first_order

  module subroutine apply_timestep_second_order(phi, diag, M, b)
    use kinds, only: ccs_int
    use mat, only: set_matrix_diagonal, get_matrix_diagonal
    use vec, only: get_vector_data, restore_vector_data
    use utils, only: update, finalise

    class(field), intent(inout) :: phi
    class(ccs_vector), intent(inout) :: diag
    class(ccs_matrix), intent(inout) :: M
    class(ccs_vector), intent(inout) :: b

    real(ccs_real), dimension(:), pointer :: diag_data
    real(ccs_real), dimension(:), pointer :: b_data
    real(ccs_real), dimension(:), pointer :: phi_old1_data
    real(ccs_real), dimension(:), pointer :: phi_old2_data
    real(ccs_real) :: rho
    integer(ccs_int) :: i
    integer(ccs_int) :: local_num_cells

    type(cell_locator) :: loc_p
    real(ccs_real) :: V_p

    rho = 1.0

    call finalise(M)
    call get_matrix_diagonal(M, diag)

    call get_vector_data(phi%old_values(1)%vec, phi_old1_data)
    call get_vector_data(phi%old_values(2)%vec, phi_old2_data)
    call get_vector_data(diag, diag_data)
    call update(b)
    call get_vector_data(b, b_data)

    call get_local_num_cells(local_num_cells)
    do i = 1, local_num_cells
      call create_cell_locator(i, loc_p)
      call get_volume(loc_p, V_p)

      ! A = A + 1.5*rho*V/dt
      diag_data(i) = diag_data(i) + 1.5 * rho * V_p / get_timestep()

      ! b = b + rho*V/dt * (2*phi_old(n-1) - 0.5*phi_old(n-2))
      b_data(i) = b_data(i) + rho * V_p / get_timestep() * (2 * phi_old1_data(i) - 0.5 * phi_old2_data(i))
    end do
    call restore_vector_data(phi%old_values(1)%vec, phi_old1_data)
    call restore_vector_data(phi%old_values(2)%vec, phi_old2_data)
    call restore_vector_data(diag, diag_data)
    call restore_vector_data(b, b_data)
    call set_matrix_diagonal(diag, M)
  end subroutine apply_timestep_second_order

  module subroutine apply_timestep_theta(theta, phi, diag, M, b)
    use kinds, only: ccs_int
    use mat, only: set_matrix_diagonal, get_matrix_diagonal
    use vec, only: get_vector_data, restore_vector_data
    use utils, only: update, finalise
    use meshing, only: get_local_num_cells

    real(ccs_real), intent(in) :: theta
    class(field), intent(inout) :: phi
    class(ccs_vector), intent(inout) :: diag
    class(ccs_matrix), intent(inout) :: M
    class(ccs_vector), intent(inout) :: b

    real(ccs_real), dimension(:), pointer :: diag_data
    real(ccs_real), dimension(:), pointer :: b_data
    real(ccs_real), dimension(:), pointer :: phi_old1_data
    real(ccs_real), dimension(:), pointer :: phi_old2_data
    real(ccs_real) :: rho
    integer(ccs_int) :: i
    integer(ccs_int) :: local_num_cells

    type(cell_locator) :: loc_p
    real(ccs_real) :: V_p

    rho = 1.0

    call finalise(M)
    call get_matrix_diagonal(M, diag)

    call get_vector_data(phi%old_values(1)%vec, phi_old1_data)
    call get_vector_data(phi%old_values(2)%vec, phi_old2_data)
    call get_vector_data(diag, diag_data)
    call update(b)
    call get_vector_data(b, b_data)

    call get_local_num_cells(local_num_cells)
    do i = 1, local_num_cells
      call create_cell_locator(i, loc_p)
      call get_volume(loc_p, V_p)

      ! A = A + (1.0 + 0.5 * theta)*rho*V/dt
      diag_data(i) = diag_data(i) + (1.0 + 0.5 * theta) * rho * V_p / get_timestep()

      ! b = b + rho*V/dt * ((1.0 + theta)*phi_old(n-1) - 0.5*theta*phi_old(n-2))
      b_data(i) = b_data(i) + rho * V_p / get_timestep() * ((1.0 + theta) * phi_old1_data(i) - 0.5 * theta * phi_old2_data(i))
    end do
    call restore_vector_data(phi%old_values(1)%vec, phi_old1_data)
    call restore_vector_data(phi%old_values(2)%vec, phi_old2_data)
    call restore_vector_data(diag, diag_data)
    call restore_vector_data(b, b_data)
    call set_matrix_diagonal(diag, M)
  end subroutine apply_timestep_theta

end submodule timestepping_common
