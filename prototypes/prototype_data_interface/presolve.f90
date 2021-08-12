module presolve
use types
implicit none

interface

  module subroutine do_initialisation(all_d)
    type(all_data), intent(inout) :: all_d
  end subroutine do_initialisation

  module subroutine initialise_turbulence(all_d)
    type(all_data), intent(inout) :: all_d
  end subroutine initialise_turbulence 

end interface

end module presolve
