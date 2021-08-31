program ccs
use stepper
use types
implicit none

type(all_data) :: all_d

call initialise_stepper(all_d)
call do_steps(all_d)

end program
