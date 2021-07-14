program ccs
use presolve
use postsolve
use flowsolve
implicit none
integer :: timestep, solvestep

call read_input()
call do_initialisation()

do timestep = 1, 3 
  call do_time_step_stuff()
  do solvestep = 1, 2
    call turbulence()
    call particles()
    call flux()
    call solve()
  end do
end do

call write_output()
call do_finalisation()
end program
