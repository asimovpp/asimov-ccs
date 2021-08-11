program ccs
use presolve
use postsolve
use flowsolve
implicit none

call read_input()
call do_initialisation()

call stepper()

call write_output()
call do_finalisation()
end program
