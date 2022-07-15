!>  Program file for scalar advection case
!
!

program bc_test
#include "ccs_macros.inc"

  !! ASiMoV-CCS uses
  use testing_lib
  use kinds, only : ccs_real, ccs_int
  use types, only : field, central_field
  use utils, only : debug_print, exit_print, str
  use boundary_conditions, only : read_bc_config

  implicit none

  class(field), allocatable :: u
  integer(ccs_int) :: i
  real(ccs_real) :: abs_err

  call init()

  ! Init velocities and scalar
  allocate(central_field :: u)

  ! Read bc configuration
  call read_bc_config("test_read_bc_config.yaml", "u", u%bcs)
  do i = 1, 4
    select case (i)
    case (1)
      if (u%bcs%names(i) /= -1) then
        call error_abort("bc name incorrect. Expected " // str(-1) // " received " // str(u%bcs%names(i)))
      end if
      if (u%bcs%bc_types(i) /= 5) then
        call error_abort("bc name incorrect. Expected " // str(5) // " received " // str(u%bcs%bc_types(i)))
      end if
      abs_err = u%bcs%values(i) - 107.345 
      if (abs_err > eps) then
        print *, u%bcs%values(i), abs_err
        call error_abort("bc name incorrect. Expected 107.345"  // " received " // str(u%bcs%values(i)))
      end if
    case (2)
      if (u%bcs%names(i) /= -4) then
        call error_abort("bc name incorrect. Expected " // str(-4) // " received " // str(u%bcs%names(i)))
      end if
      if (u%bcs%bc_types(i) /= 1) then
        call error_abort("bc name incorrect. Expected " // str(1) // " received " // str(u%bcs%bc_types(i)))
      end if
    case (3)
      if (u%bcs%names(i) /= -2) then
        call error_abort("bc name incorrect. Expected " // str(-2) // " received " // str(u%bcs%names(i)))
      end if
      if (u%bcs%bc_types(i) /= 3) then
        call error_abort("bc name incorrect. Expected " // str(3) // " received " // str(u%bcs%bc_types(i)))
      end if
      abs_err = abs(u%bcs%values(i) - (-5))
      if (abs_err > eps) then
        call error_abort("bc name incorrect. Expected " // str(-5) // " received " // str(u%bcs%values(i)))
      end if
    case (4)
      if (u%bcs%names(i) /= -3) then
        call error_abort("bc name incorrect. Expected " // str(-3) // " received " // str(u%bcs%names(i)))
      end if
      if (u%bcs%bc_types(i) /= 1) then
        call error_abort("bc name incorrect. Expected " // str(1) // " received " // str(u%bcs%bc_types(i)))
      end if
    end select
  end do

  call fin()
end program bc_test
