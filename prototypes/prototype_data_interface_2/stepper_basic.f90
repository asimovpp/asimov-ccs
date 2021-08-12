submodule (stepper) stepper_2loop
use problem_former
use solver
implicit none

contains
  
  module subroutine initialise_stepper(d)
    type(all_data), intent(inout) :: d
    call initialise_form_problem(d)
  end subroutine initialise_stepper
  
  module subroutine do_steps(d)
    integer :: timestep, solvestep
    type(all_data), intent(inout) :: d

    print *,"Using basic stepper"
    do timestep = 1, 3 
      do solvestep = 1, 2
        call form_problem(d%form_problem_d)
        call solve()
      end do
    end do
  end subroutine do_steps
end submodule stepper_2loop

