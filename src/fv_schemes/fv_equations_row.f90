!v Submodule file fv_equations_row.smod
submodule(fv_equations) fv_equations_row

  use utils, only: set_col, set_entry, set_row

  implicit none

contains

  module subroutine equation_prepare_row(self, global_row_index, required_capacity)
    class(equation), intent(inout) :: self
    integer(ccs_int), intent(in) :: global_row_index
    integer(ccs_int), intent(in), optional :: required_capacity

    integer(ccs_int) :: min_capacity

    min_capacity = 0_ccs_int
    if (present(required_capacity)) min_capacity = required_capacity

    if (global_row_index <= 0_ccs_int) then
      error stop "equation%prepare_row: invalid global row index"
    end if

    call ensure_row_capacity(self%row, min_capacity)

    self%row%global_row_index = global_row_index
    self%row%n_entries = 0_ccs_int
    self%row%rhs = 0.0_ccs_real
    self%row%is_ready = .true.

  end subroutine equation_prepare_row

  module subroutine equation_add_row_entry(self, global_col_index, coefficient)
    class(equation), intent(inout) :: self
    integer(ccs_int), intent(in) :: global_col_index
    real(ccs_real), intent(in) :: coefficient

    integer(ccs_int) :: entry

    if (.not. self%row%is_ready) then
      error stop "equation%add_row_entry: row not prepared"
    end if

    entry = self%row%n_entries + 1_ccs_int
    call ensure_row_capacity(self%row, entry)

    self%row%n_entries = entry
    self%row%global_col_indices(entry) = global_col_index
    self%row%coefficients(entry) = coefficient

  end subroutine equation_add_row_entry

  module subroutine equation_set_rhs(self, rhs)
    class(equation), intent(inout) :: self
    real(ccs_real), intent(in) :: rhs

    if (.not. self%row%is_ready) then
      error stop "equation%set_rhs: row not prepared"
    end if

    self%row%rhs = rhs

  end subroutine equation_set_rhs

  module subroutine equation_flush_row(self, mat_coeffs, vec_values)
    class(equation), intent(in) :: self
    type(matrix_values), intent(inout) :: mat_coeffs
    type(vector_values), intent(inout) :: vec_values

    integer(ccs_int) :: i

    call check_flushable_row(self%row, "equation%flush_row")

    call set_row(self%row%global_row_index, mat_coeffs)
    do i = 1, self%row%n_entries
      call set_col(self%row%global_col_indices(i), mat_coeffs)
      call set_entry(self%row%coefficients(i), mat_coeffs)
    end do

    call self%flush_rhs(vec_values)

  end subroutine equation_flush_row

  module subroutine equation_flush_rhs(self, vec_values)
    class(equation), intent(in) :: self
    type(vector_values), intent(inout) :: vec_values

    call check_flushable_row(self%row, "equation%flush_rhs")

    call set_row(self%row%global_row_index, vec_values)
    call set_entry(self%row%rhs, vec_values)

  end subroutine equation_flush_rhs

  module subroutine ensure_row_capacity(row, required)
    type(equation_row), intent(inout) :: row
    integer(ccs_int), intent(in) :: required

    integer(ccs_int) :: capacity
    integer(ccs_int) :: extension_size
    integer(ccs_int), allocatable :: global_col_indices_extension(:)
    real(ccs_real), allocatable :: coefficients_extension(:)

    if (required < 0_ccs_int) then
      error stop "ensure_row_capacity: invalid required capacity"
    end if
    if (required == 0_ccs_int) return

    if (.not. allocated(row%global_col_indices) .and. .not. allocated(row%coefficients)) then
      allocate (row%global_col_indices(required), source=0_ccs_int)
      allocate (row%coefficients(required), source=0.0_ccs_real)
      return
    end if
    if (.not. allocated(row%global_col_indices) .or. .not. allocated(row%coefficients)) then
      error stop "ensure_row_capacity: inconsistent row storage"
    end if

    capacity = size(row%global_col_indices, kind=ccs_int)
    if (size(row%coefficients, kind=ccs_int) /= capacity) then
      error stop "ensure_row_capacity: inconsistent row storage"
    end if
    if (capacity >= required) return

    extension_size = required - capacity
    allocate (global_col_indices_extension(extension_size), source=0_ccs_int)
    allocate (coefficients_extension(extension_size), source=0.0_ccs_real)
    row%global_col_indices = [row%global_col_indices, global_col_indices_extension]
    row%coefficients = [row%coefficients, coefficients_extension]

  end subroutine ensure_row_capacity

  subroutine check_flushable_row(row, caller)
    type(equation_row), intent(in) :: row
    character(len=*), intent(in) :: caller

    if (.not. row%is_ready) then
      error stop caller // ": row not prepared"
    end if
    if (row%global_row_index <= 0_ccs_int) then
      error stop caller // ": invalid global row index"
    end if
    if (row%n_entries < 0_ccs_int) then
      error stop caller // ": invalid entry count"
    end if
    if (row%n_entries == 0_ccs_int) return

    if (.not. allocated(row%global_col_indices) .or. .not. allocated(row%coefficients)) then
      error stop caller // ": row entries not allocated"
    end if
    if (size(row%global_col_indices, kind=ccs_int) < row%n_entries .or. &
        size(row%coefficients, kind=ccs_int) < row%n_entries) then
      error stop caller // ": row entry storage too small"
    end if

  end subroutine check_flushable_row

end submodule fv_equations_row
