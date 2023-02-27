module custom_matrix

  use kinds, only: ccs_int, ccs_real

  implicit none

  !v a "csr-lite" matrix format that assumes a constant number
  !  of non-zeros in each row
  type, public :: csr_matrix
    real(ccs_real), dimension(:), allocatable :: values
    integer(ccs_int), dimension(:), allocatable :: columns
    integer(ccs_int) :: values_per_row
  end type csr_matrix

contains

  subroutine create_new_matrix(rows, values_per_row, matrix)

    integer(ccs_int), intent(in) :: rows
    integer(ccs_int), intent(in) :: values_per_row
    type(csr_matrix), intent(inout) :: matrix

    matrix%values_per_row = values_per_row

    allocate(matrix%values(rows * values_per_row))
    matrix%values(:) = 0.0

    allocate(matrix%columns(rows * values_per_row))
    matrix%columns(:) = 0.0

  end subroutine

  subroutine insert_values(row, vals, cols, matrix)

    integer(ccs_int), intent(in) :: row
    real(ccs_real), dimension(:) :: vals
    integer(ccs_int), dimension(:) :: cols
    type(csr_matrix), intent(inout) :: matrix

    integer(ccs_int) :: offset
    
    offset = (row - 1) * matrix%values_per_row + 1

    matrix%values(offset : offset + matrix%values_per_row) = vals(:)
    matrix%columns(offset : offset + matrix%values_per_row) = cols(:)

  end subroutine

  subroutine print_matrix(matrix)

    type(csr_matrix), intent(in) :: matrix

    print *, "values_per_row = ", matrix%values_per_row
    print *, "values"
    print *, matrix%values
    print *, "columns"
    print *, matrix%columns

  end subroutine

end module custom_matrix




program test_matrix
  use custom_matrix

  implicit none

  type(csr_matrix) :: m

  call create_new_matrix(3, 4, m)
  call print_matrix(m)
  call insert_values(1, (/1.0_ccs_real,2.0_ccs_real,3.0_ccs_real,4.0_ccs_real/), (/1,2,3,4/), m)
  call insert_values(2, (/5.0_ccs_real,6.0_ccs_real,7.0_ccs_real,8.0_ccs_real/), (/2,3,4,1/), m)
  call insert_values(3, (/9.0_ccs_real,10.0_ccs_real,11.0_ccs_real,12.0_ccs_real/), (/3,4,1,2/), m)
  call print_matrix(m)

end program
