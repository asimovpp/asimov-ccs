module my_math
implicit none
real, parameter :: pi = 4.*atan(1.)
real :: tau

interface
  module subroutine pi_mult(pi,tau)
    real, intent(in) :: pi
    real, intent(out) :: tau
  end subroutine pi_mult
end interface
contains
end module my_math
