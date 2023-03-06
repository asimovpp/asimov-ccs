module custom_vector

  use kinds, only: ccs_int, ccs_real

  implicit none

  !v a simple vector format
  type, public :: custom_vector
    real(ccs_real), dimension(:), allocatable :: values
    integer(ccs_int) :: num_rows
  end type custom_vector

contains

  subroutine create_new_vector(rows, vector)

    integer(ccs_int), intent(in) :: rows
    type(custom_vector), intent(inout) :: vector

    vector%num_rows = rows

    allocate(vector%values(rows))
    vector%values(:) = 0.0


  end subroutine

  subroutine zeros(vector)

    type(custom_vector), intent(inout) :: vector
    
    vector%values(:) = 0.0_ccs_real

  end subroutine

  subroutine add_values(row, val, vector)

    integer(ccs_int), intent(in) :: row
    real(ccs_real) :: val
    type(custom_vector), intent(inout) :: vector

    vector%values(row) = vector%values(row) + val

  end subroutine

  subroutine insert_values(row, val, vector)

    integer(ccs_int), intent(in) :: row
    real(ccs_real) :: val
    type(custom_vector), intent(inout) :: vector

    vector%values(row) = val
  end subroutine

  subroutine print_vector(vector)

    type(custom_vector), intent(in) :: vector

    print *, "values"
    print *, vector%values
    print *, "num row"
    print *, vector%num_rows

  end subroutine

end module custom_vector




! program test_matrix
!   use custom_matrix
! 
!   implicit none
! 
!   type(csr_matrix) :: m
! 
!   call create_new_matrix(3, 4, m)
!   call print_matrix(m)
!   call insert_values(1, (/1.0_ccs_real,2.0_ccs_real,3.0_ccs_real,4.0_ccs_real/), (/1,2,3,4/), m)
!   call insert_values(2, (/5.0_ccs_real,6.0_ccs_real,7.0_ccs_real,8.0_ccs_real/), (/2,3,4,1/), m)
!   call insert_values(3, (/9.0_ccs_real,10.0_ccs_real,11.0_ccs_real,12.0_ccs_real/), (/3,4,1,2/), m)
!   call print_matrix(m)
! 
! end program
! 