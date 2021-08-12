module solver_types
use type_headers
implicit none

type, extends(solver_data) :: solver_amg_data
  real :: x1
  real, pointer :: datar(:)
  !real, allocatable :: datar(:)
end type solver_amg_data 

type, extends(solver_data) :: solver_cgstab_data
  real :: x1
end type solver_cgstab_data 

end module
