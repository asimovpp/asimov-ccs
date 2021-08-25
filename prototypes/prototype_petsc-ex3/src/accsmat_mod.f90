module accsmat

  use accs_types, only : matrix, matrix_init_data, matrix_values
  
  implicit none

  private

  public :: create_matrix, update_matrix, begin_update_matrix, end_update_matrix, set_matrix_values

  interface
     module subroutine create_matrix(mat_dat, M)
       type(matrix_init_data), intent(in) :: mat_dat
       class(matrix), allocatable, intent(out) :: M
     end subroutine

     module subroutine update_matrix(M)
       class(matrix), intent(inout) :: M
     end subroutine
     module subroutine begin_update_matrix(M)
       class(matrix), intent(inout) :: M
     end subroutine
     module subroutine end_update_matrix(M)
       class(matrix), intent(inout) :: M
     end subroutine

     module subroutine set_matrix_values(mat_values, M)
       type(matrix_values), intent(in) :: mat_values
       class(matrix), intent(inout) :: M
     end subroutine
  end interface
  
end module accsmat
