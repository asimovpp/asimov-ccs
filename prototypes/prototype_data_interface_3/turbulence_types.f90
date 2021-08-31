module turbulence_types
use type_headers
implicit none

type, extends(turbulence_data) :: turbulence_data_1
  real :: c1
end type turbulence_data_1
type, extends(turbulence_data) :: turbulence_data_2
  real :: c2
  real :: c3
end type turbulence_data_2

end module
