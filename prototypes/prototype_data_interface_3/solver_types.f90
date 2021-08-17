module solver_types
use type_headers
use problem_former_types
implicit none

type, extends(solver_data) :: solver_amg_data
  real :: x1
  type(form_problem_data_basic), pointer :: fpd
end type solver_amg_data 

type, extends(solver_data) :: solver_cgstab_data
  real :: x1
end type solver_cgstab_data 

end module
