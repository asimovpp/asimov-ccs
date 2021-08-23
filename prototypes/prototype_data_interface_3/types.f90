module types
use type_headers
use turbulence_types
use problem_former_types
use solver_types
implicit none

type :: all_data
  class(turbulence_data), allocatable :: turbulence_d
  class(form_problem_data), allocatable :: form_problem_d
  class(solver_data), allocatable :: solver_d
end type

end module
