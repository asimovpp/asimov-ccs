module custom_vector

  use kinds, only: ccs_int, ccs_real

  implicit none

  !v a simple vector format
  type, public :: lite_vector
    real(ccs_real), dimension(:), allocatable :: values
    integer(ccs_int) :: num_rows
  end type lite_vector

contains

  subroutine create_new_vector(rows, vector)

    integer(ccs_int), intent(in) :: rows
    type(lite_vector), intent(inout) :: vector

    vector%num_rows = rows

    allocate(vector%values(rows))
    vector%values(:) = 0.0

  end subroutine

  subroutine zeros(vector)

    type(lite_vector), intent(inout) :: vector
    
    vector%values(:) = 0.0_ccs_real

  end subroutine

  subroutine add_value(row, val, vector)

    integer(ccs_int), intent(in) :: row
    real(ccs_real), intent(in) :: val
    type(lite_vector), intent(inout) :: vector

    vector%values(row) = vector%values(row) + val

  end subroutine

  subroutine insert_value(row, val, vector)

    integer(ccs_int), intent(in) :: row
    real(ccs_real), intent(in):: val
    type(lite_vector), intent(inout) :: vector

    vector%values(row) = val
  end subroutine

  subroutine print_vector(vector)

    type(lite_vector), intent(in) :: vector

    print *, "values"
    print *, vector%values
    print *, "num row"
    print *, vector%num_rows

  end subroutine

end module custom_vector




program test_vector
  use custom_vector

  implicit none

  type(lite_vector) :: v

  call create_new_vector(4, v)
  call print_vector(v)
  call insert_value(1, 1.0_ccs_real, v)
  call insert_value(2, 1.0_ccs_real, v)
  call insert_value(3, 3.0_ccs_real, v)
  call insert_value(4, 4.0_ccs_real, v)
  call print_vector(v)

  call add_value(2, 1.0_ccs_real, v)
  call print_vector(v)

  call zeros(v)
  call print_vector(v)

end program
