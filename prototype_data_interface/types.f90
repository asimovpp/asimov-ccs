module types
implicit none

type :: turbulence_input
end type turbulence_input
type, extends(turbulence_input) :: turb_1
  real :: c1
end type turb_1
type, extends(turbulence_input) :: turb_2
  real :: c2
  real :: c3
end type turb_2

type :: all_data
  class(turbulence_input), allocatable :: turb_data
end type

end module
