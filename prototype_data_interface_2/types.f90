module types
implicit none

type :: turbulence_data
end type turbulence_data
type, extends(turbulence_data) :: turbulence_data_1
  real :: c1
end type turbulence_data_1
type, extends(turbulence_data) :: turbulence_data_2
  real :: c2
  real :: c3
end type turbulence_data_2

type :: form_problem_data
end type form_problem_data
type, extends(form_problem_data) :: form_problem_data_basic
  real :: x1
end type form_problem_data_basic

type :: all_data
  class(turbulence_data), allocatable :: turbulence_d
  class(form_problem_data), allocatable :: form_problem_d
end type

end module
