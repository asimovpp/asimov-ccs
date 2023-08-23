program test_bc_profile_interpolation
#include "ccs_macros.inc"

  use testing_lib
  use kinds, only: ccs_real, ccs_int
  use fv, only: get_value_from_bc_profile

  implicit none

  type(bc_profile), allocatable :: profile
  real(ccs_real), dimension(3) :: x
  real(ccs_real) :: bc_value, expected_value
  integer(ccs_int) :: i


  x = (/ 0.0_ccs_real, 0.0_ccs_real, 0.0_ccs_real/)

  call init()

  call generate_profile(profile)
  print *, " Coordinates: ", profile%coordinates
  print *, " Values     : ", profile%values

  x(1) = 0.0_ccs_real
  call get_value_from_bc_profile(x, profile, bc_value)

  call assert_eq(bc_value, 0.0_ccs_real, "wrong lower bound interpolation value")

  do i=1, 5
    x(1) = real(i, ccs_real) * 1.5_ccs_real + 0.2_ccs_real
    expected_value = f(x(1))
    call get_value_from_bc_profile(x, profile, bc_value)

    call assert_eq(bc_value, expected_value, "wrong in bound interpolation value")
  end do

  x(1) = 12.0_ccs_real
  call get_value_from_bc_profile(x, profile, bc_value)

  call assert_eq(bc_value, 20.0_ccs_real, "wrong higher bound interpolation value")

  call fin()

contains

function f(x) result(y)
  real(ccs_real), intent(in) :: x
  real(ccs_real) :: y

  if (x .lt. 5.0_ccs_real) then
    y = x
  else
    y = 5.0_ccs_real + 3_ccs_real * (x - 5.0_ccs_real)
  end if

end function

subroutine generate_profile(profile)

  type(bc_profile), allocatable, intent(out) :: profile
  integer(ccs_int) :: i

  allocate(profile)
  allocate(profile%centre(3))
  allocate(profile%values(11))
  allocate(profile%coordinates(11))

  profile%centre = (/ 0.0_ccs_real, 0.0_ccs_real, 0.0_ccs_real/)

  do i=0, 10
    profile%coordinates(i+1) = real(i, ccs_real)
    profile%values(i+1) = f(real(i, ccs_real))
  end do

end subroutine

end program test_bc_profile_interpolation