module timestepping_common_types

  use kinds, only: ccs_real, ccs_int

  implicit none

  private
  public :: ptr_handle

  ! Small type used to enable creating an array of data pointers, this allows the timestep
  ! application to be implemented more generically across arbitrary transient stencil
  ! widths. Additionally this allows for data hiding.
  type ptr_handle
     private
     real(ccs_real), dimension(:), pointer :: data => null()
   contains
     procedure :: read => read_ptr_handle
     procedure :: get_pointer
     procedure :: set_pointer
  end type ptr_handle

contains
  
  pure real(ccs_real) function read_ptr_handle(self, idx) result(val)
    class(ptr_handle), intent(in) :: self
    integer(ccs_int), intent(in) :: idx

    val = self%data(idx)
    
  end function read_ptr_handle

  function get_pointer(self) result(ptr)
    class(ptr_handle), intent(in) :: self
    real(ccs_real), dimension(:), pointer :: ptr

    ptr => self%data
  end function get_pointer

  subroutine set_pointer(self, ptr)
    class(ptr_handle), intent(out) :: self
    real(ccs_real), dimension(:), pointer :: ptr

    self%data => ptr
  end subroutine set_pointer
  
end module timestepping_common_types

submodule(timestepping) timestepping_common
#include "ccs_macros.inc"

  use meshing, only: get_local_num_cells, create_cell_locator, get_volume
  use transient_kernels, only: transient_kernel
  use types, only: cell_locator
  use utils, only: exit_print

  use timestepping_common_types, only: ptr_handle

  implicit none

  logical :: timestepping_active = .false. !< flag to signify whether timestepping should occur
  logical :: timestep_is_set = .false. !< flag to signify whether dt has already been set
  real(ccs_real) :: dt = huge(0.0_ccs_real) !< timestep size
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

    integer(ccs_int), intent(in) :: num_old_vals
    class(field), intent(inout) :: x

    integer(ccs_int) :: i

    do i = num_old_vals, 2, -1
      call copy_old_data(x%old_values(i - 1)%vec, x%old_values(i)%vec)
    end do
    call copy_old_data(x%values, x%old_values(1)%vec)

  end subroutine

  !> Copies vector data from new to old destination using readonly views of current data.
  subroutine copy_old_data(newv, oldv)
    use vec, only: get_vector_data, restore_vector_data, get_vector_data_readonly, restore_vector_data_readonly

    class(ccs_vector), intent(inout) :: newv
    class(ccs_vector), intent(inout) :: oldv

    real(ccs_real), dimension(:), pointer :: new_data, old_data

    call get_vector_data_readonly(newv, new_data)
    call get_vector_data(oldv, old_data)

    old_data(:) = new_data(:)

    call restore_vector_data_readonly(newv, new_data)
    call restore_vector_data(oldv, old_data)
    
  end subroutine copy_old_data

  module subroutine apply_timestep_kernel(transient, phi, diag, M, b)
    use kinds, only: ccs_int
    use mat, only: set_matrix_diagonal, get_matrix_diagonal
    use vec, only: get_vector_data, restore_vector_data, get_vector_data_readonly, restore_vector_data_readonly
    use utils, only: update, finalise

    class(transient_kernel), intent(inout) :: transient ! The transient kernel
    class(field), intent(inout) :: phi
    class(ccs_vector), intent(inout) :: diag
    class(ccs_matrix), intent(inout) :: M
    class(ccs_vector), intent(inout) :: b

    real(ccs_real), dimension(:), pointer :: diag_data
    real(ccs_real), dimension(:), pointer :: b_data
    real(ccs_real), dimension(:), pointer :: ptr
    type(ptr_handle), dimension(:), allocatable :: old_pointer
    integer(ccs_int) :: i

    call transient%init()
    call transient%set_step(current_step+1)
    call transient%set_dt(get_timestep())

    call finalise(M)
    call get_matrix_diagonal(M, diag)

    allocate(old_pointer(transient%get_width()))
    do i = 1, transient%get_width()
       call get_vector_data_readonly(phi%old_values(i)%vec, ptr)
       call old_pointer(i)%set_pointer(ptr)
       nullify(ptr)
    end do

    call get_vector_data(diag, diag_data)
    call update(b)
    call get_vector_data(b, b_data)

    call apply_kernel_driver(transient, old_pointer, diag_data, b_data)

    do i = 1, transient%get_width()
       ptr => old_pointer(i)%get_pointer()
       call restore_vector_data_readonly(phi%old_values(i)%vec, ptr)
       nullify(ptr)
    end do
    call restore_vector_data(diag, diag_data)
    call restore_vector_data(b, b_data)
    call set_matrix_diagonal(diag, M)

  end subroutine apply_timestep_kernel

  subroutine apply_kernel_driver(transient, old_pointer, diag_data, b_data)

    class(transient_kernel), intent(inout) :: transient
    type(ptr_handle), dimension(:), intent(in) :: old_pointer
    real(ccs_real), dimension(:), intent(inout) :: diag_data, b_data

    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: i, j
    type(cell_locator) :: loc_p
    real(ccs_real) :: V_p, coeff, rhs
    real(ccs_real) :: rho

    real(ccs_real), allocatable, dimension(:) :: old

    rho = 1.0

    allocate(old(transient%get_width()))
    
    call get_local_num_cells(local_num_cells)
    !$omp parallel do default(none) schedule(static)   &
    !$omp shared(local_num_cells, old_pointer, rho,    &
    !$omp transient, diag_data, b_data)                &
    !$omp private(i, j, old, loc_p, V_p, coeff, rhs)
    do i = 1, local_num_cells
      call create_cell_locator(i, loc_p)
      call get_volume(loc_p, V_p)

      call transient%eval_coeffs(rho, V_p, coeff)

      do j = 1, transient%get_width()
         old(j) = old_pointer(j)%read(i)
      end do

      call transient%eval_explicit(rho, V_p, old, rhs)

      diag_data(i) = diag_data(i) + coeff
      b_data(i) = b_data(i) + rhs
    end do
    !$omp end parallel do

  end subroutine apply_kernel_driver

end submodule timestepping_common
