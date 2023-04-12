submodule(timestepping) timestepping_common
#include "ccs_macros.inc"

  use utils, only: exit_print

  implicit none

  logical :: timestepping_active = .false. !< flag to signify whether timestepping should occur
  logical :: timestep_is_set = .false. !< flag to signify whether dt has already been set
  real(ccs_real) :: dt !< timestep size

contains
  
  module subroutine activate_timestepping()
    timestepping_active = .true.
  end subroutine

  module subroutine apply_timestep_first_order(mesh, phi, diag, M, b)

    use kinds, only: ccs_int
    use mat, only: set_matrix_diagonal, get_matrix_diagonal
    use vec, only: get_vector_data, restore_vector_data
    use utils, only: update, finalise
    use meshing, only: get_local_num_cells

    type(ccs_mesh), intent(in) :: mesh
    class(field), intent(inout) :: phi
    class(ccs_vector), intent(inout) :: diag
    class(ccs_matrix), intent(inout) :: M
    class(ccs_vector), intent(inout) :: b

    real(ccs_real), dimension(:), pointer :: diag_data
    real(ccs_real), dimension(:), pointer :: b_data
    real(ccs_real), dimension(:), pointer :: phi_data
    integer(ccs_int) :: i
    integer(ccs_int) :: local_num_cells

    ! V = mesh%volumes
    call finalise(M)
    call get_matrix_diagonal(M, diag)

    call get_vector_data(phi%old_values(1)%vec, phi_data)
    call get_vector_data(diag, diag_data)
    call update(b)
    call get_vector_data(b, b_data)

    call get_local_num_cells(mesh, local_num_cells)
    do i = 1, local_num_cells
      ! A = A + V/dt
      diag_data(i) = diag_data(i) + mesh%geo%volumes(i) / get_timestep()

      ! b = b + V/dt * phi_old
      b_data(i) = b_data(i) + mesh%geo%volumes(i) / get_timestep() * phi_data(i)
    end do
    call restore_vector_data(phi%old_values(1)%vec, phi_data)
    call restore_vector_data(diag, diag_data)
    call restore_vector_data(b, b_data)
    call set_matrix_diagonal(diag, M)

  end subroutine apply_timestep_first_order
  
  module subroutine apply_timestep_second_order(mesh, phi, diag, M, b)
    use kinds, only: ccs_int
    use mat, only: set_matrix_diagonal, get_matrix_diagonal
    use vec, only: get_vector_data, restore_vector_data
    use utils, only: update, finalise
    use meshing, only: get_local_num_cells

    type(ccs_mesh), intent(in) :: mesh
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

    rho = 1.0

    ! V = mesh%volumes
    call finalise(M)
    call get_matrix_diagonal(M, diag)

    call get_vector_data(phi%old_values(1)%vec, phi_old1_data)
    call get_vector_data(phi%old_values(2)%vec, phi_old2_data)
    call get_vector_data(diag, diag_data)
    call update(b)
    call get_vector_data(b, b_data)

    call get_local_num_cells(mesh, local_num_cells)
    do i = 1, local_num_cells
      ! A = A + 1.5*rho*V/dt
      diag_data(i) = diag_data(i) + 1.5 * rho * mesh%geo%volumes(i) / get_timestep()

      ! b = b + rho*V/dt * (2*phi_old(n-1) - 0.5*phi_old(n-2))
      b_data(i) = b_data(i) + rho * mesh%geo%volumes(i) / get_timestep() * (2 * phi_old1_data(i) - 0.5 * phi_old2_data(i))
    end do
    call restore_vector_data(phi%old_values(1)%vec, phi_old1_data)
    call restore_vector_data(phi%old_values(2)%vec, phi_old2_data)
    call restore_vector_data(diag, diag_data)
    call restore_vector_data(b, b_data)
    call set_matrix_diagonal(diag, M)
  end subroutine apply_timestep_second_order

  module function timestepping_is_active() result(active)
    logical :: active
    active = timestepping_active
  end function
  
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
  
  module subroutine update_old_values_generic(num_old_vals, x)

    use vec, only: get_vector_data, restore_vector_data

    integer(ccs_int), intent(in) :: num_old_vals
    class(field), intent(inout) :: x

    real(ccs_real), dimension(:), pointer :: values_data, old_values_data
    integer(ccs_int) :: i
    
    do i = num_old_vals, 2, -1
      call get_vector_data(x%old_values(i)%vec, old_values_data)
      call get_vector_data(x%old_values(i-1)%vec, values_data)
      old_values_data = values_data
      call restore_vector_data(x%old_values(i)%vec, old_values_data)
      call restore_vector_data(x%old_values(i-1)%vec, values_data)
    end do

    call get_vector_data(x%old_values(1)%vec, old_values_data)
    call get_vector_data(x%values, values_data)
    old_values_data = values_data
    call restore_vector_data(x%old_values(1)%vec, old_values_data)
    call restore_vector_data(x%values, values_data)

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


end submodule timestepping_common
