module solver
use types
implicit none

interface
  
  module subroutine initialise_solve(d)
    type(all_data), intent(inout) :: d
  end subroutine initialise_solve

  module subroutine solve(d)
    class(solver_data), intent(inout) :: d
  end subroutine solve

end interface

end module solver
