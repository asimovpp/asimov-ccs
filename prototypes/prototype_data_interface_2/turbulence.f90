module turbulence
use types
implicit none

interface
  module subroutine initialise_turbulence(d)
    type(all_data), intent(inout) :: d
  end subroutine initialise_turbulence 
  
  module subroutine do_turbulence(d)
    class(turbulence_data), intent(inout) :: d
  end subroutine do_turbulence 

end interface

end module turbulence
