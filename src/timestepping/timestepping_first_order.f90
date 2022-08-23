submodule (timestepping) timestepping_first_order
#include "ccs_macros.inc"

  use utils, only: exit_print
  
  implicit none

  real(ccs_real) :: dt
  logical :: timestep_is_set = .false.
  logical :: timestepping_is_active = .false.

contains

  module subroutine apply_timestep(mesh, phi, diag, M, b)
    use kinds, only: ccs_int
    use mat, only : set_matrix_diagonal, get_matrix_diagonal
    use vec, only: get_vector_data, restore_vector_data
    use utils, only: update, finalise
    
    type(ccs_mesh), intent(in) :: mesh
    class(field), intent(in) :: phi
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
   
    call get_vector_data(phi%old_values, phi_data)
    call get_vector_data(diag, diag_data)
    call update(b)
    call get_vector_data(b, b_data)

    
    do i = 1, mesh%nlocal
    ! A = A + V/dt
      diag_data(i) = diag_data(i) + mesh%volumes(i) / get_timestep()

    ! b = b + V/dt * phi_old
      b_data(i) = b_data(i) + mesh%volumes(i) / get_timestep() * phi_data(i)
    end do
    call restore_vector_data(phi%old_values, phi_data)
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

    !u%old_values = x%values
    call get_vector_data(x%values, values_data)
    call get_vector_data(x%old_values, old_values_data)
    old_values_data = values_data
    call restore_vector_data(x%old_values, old_values_data)
    call restore_vector_data(x%values, values_data)
  end subroutine

  module subroutine activate_timestepping()
    timestepping_is_active = .true.
  end subroutine

end submodule 
