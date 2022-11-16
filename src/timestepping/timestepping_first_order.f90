!v Implementation of a first order timestepping scheme with a fixed-size timestep.
submodule(timestepping) timestepping_first_order
#include "ccs_macros.inc"

  use utils, only: exit_print

  implicit none

  real(ccs_real) :: dt !< timestep size
  logical :: timestep_is_set = .false. !< flag to signify whether dt has already been set
  logical :: timestepping_is_active = .false. !< flag to signify whether timestepping should occur
  integer(ccs_int), parameter :: num_old_vals = 1 !< the number of old field values the scheme uses

contains

  module subroutine finalise_timestep()
    !!! NOTHING TO SEE HERE!
  end subroutine finalise_timestep

  module subroutine apply_timestep(mesh, phi, diag, M, b)
    use kinds, only: ccs_int
    use mat, only: set_matrix_diagonal, get_matrix_diagonal
    use vec, only: get_vector_data, restore_vector_data
    use utils, only: update, finalise

    type(ccs_mesh), intent(in) :: mesh
    class(field), intent(inout) :: phi
    class(ccs_vector), intent(inout) :: diag
    class(ccs_matrix), intent(inout) :: M
    class(ccs_vector), intent(inout) :: b

    real(ccs_real), dimension(:), pointer :: diag_data
    real(ccs_real), dimension(:), pointer :: b_data
    real(ccs_real), dimension(:), pointer :: phi_data
    integer(ccs_int) :: i

    if (.not. timestepping_is_active) then
      return
    end if

    ! V = mesh%volumes
    call finalise(M)
    call get_matrix_diagonal(M, diag)

    call get_vector_data(phi%old_values(1)%vec, phi_data)
    call get_vector_data(diag, diag_data)
    call update(b)
    call get_vector_data(b, b_data)

    do i = 1, mesh%topo%local_num_cells
      ! A = A + V/dt
      diag_data(i) = diag_data(i) + mesh%geo%volumes(i) / get_timestep()

      ! b = b + V/dt * phi_old
      b_data(i) = b_data(i) + mesh%geo%volumes(i) / get_timestep() * phi_data(i)
    end do
    call restore_vector_data(phi%old_values(1)%vec, phi_data)
    call restore_vector_data(diag, diag_data)
    call restore_vector_data(b, b_data)
    call set_matrix_diagonal(diag, M)

  end subroutine apply_timestep

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

  module subroutine update_old_values(x)
    use vec, only: get_vector_data, restore_vector_data

    class(field), intent(inout) :: x

    real(ccs_real), dimension(:), pointer :: values_data, old_values_data

    if (timestepping_is_active) then
      !u%old_values(1) = x%values
      call get_vector_data(x%values, values_data)
      call get_vector_data(x%old_values(1)%vec, old_values_data)
      old_values_data = values_data
      call restore_vector_data(x%old_values(1)%vec, old_values_data)
      call restore_vector_data(x%values, values_data)
    end if

  end subroutine

  module subroutine activate_timestepping()
    timestepping_is_active = .true.
  end subroutine

  module subroutine initialise_old_values(vec_properties, x)
    use types, only: vector_spec
    use vec, only: create_vector
    use utils, only: update

    type(vector_spec), intent(in) :: vec_properties
    class(field), intent(inout) :: x

    if (.not. allocated(x%old_values)) then
      allocate (x%old_values(num_old_vals))
    end if

    call create_vector(vec_properties, x%old_values(1)%vec)
    call update(x%old_values(1)%vec)

  end subroutine

end submodule
