module problem_former_types
use type_headers
implicit none

type, extends(form_problem_data) :: form_problem_data_basic
  real :: x1
  real, allocatable :: datar(:)
contains
  procedure :: get_data => get_datar
end type form_problem_data_basic

contains
  real function get_datar(this,i)
    implicit none
    class(form_problem_data_basic), intent(in) :: this
    integer, intent(in) :: i
    get_datar = this%datar(i)
  end function get_datar


end module
