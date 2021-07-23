program ccs
use presolve
use flowsolve
use types
implicit none

type(all_data) :: all_d

call do_initialisation(all_d)
call stepper(all_d)

end program
