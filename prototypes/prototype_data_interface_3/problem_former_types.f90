module problem_former_types
use type_headers
implicit none

type, extends(form_problem_data) :: form_problem_data_basic
  real :: x1
  real, allocatable :: datar(:)
end type form_problem_data_basic

end module
