module flowsolve
use types
implicit none

interface
  module subroutine stepper(all_d)
    type(all_data), intent(inout) :: all_d
  end subroutine stepper 

  module subroutine turbulence(input)
    class(turbulence_input), intent(inout) :: input
  end subroutine turbulence 
  
  module subroutine solve()
  end subroutine solve
end interface

end module flowsolve
