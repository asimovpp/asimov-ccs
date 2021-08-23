module accsmat

  use accs_types, only : matrix, matrix_init_data
  
  implicit none

  private

  public :: create_matrix

  interface
     module subroutine create_matrix(mat_dat, M)
       type(matrix_init_data), intent(in) :: mat_dat
       class(matrix), allocatable, intent(out) :: M
     end subroutine
  end interface
  
end module accsmat
