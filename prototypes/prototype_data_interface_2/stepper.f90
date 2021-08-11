module stepper
use types
implicit none

interface
  
  module subroutine initialise_stepper(d)
    type(all_data), intent(inout) :: d
  end subroutine initialise_stepper

  module subroutine do_steps(d)
    type(all_data), intent(inout) :: d
  end subroutine do_steps

end interface

end module stepper
