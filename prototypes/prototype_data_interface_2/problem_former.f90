module problem_former
use types
implicit none

interface
  
  module subroutine initialise_form_problem(d)
    type(all_data), intent(inout) :: d
  end subroutine initialise_form_problem

  module subroutine form_problem(d)
    class(form_problem_data), intent(inout) :: d
  end subroutine form_problem
  
end interface

end module problem_former
