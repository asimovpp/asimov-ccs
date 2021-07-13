program ccs
use presolve
use postsolve
use flowsolve
implicit none
integer :: continue_timestep, continue_solve

call read_input()
call do_initialisation()
continue_timestep = 1
continue_solve = 1

do while (continue_timestep .eq. 1)
  call do_time_step_stuff()
  do while (continue_solve .eq. 1)
    call turbulence()
    call particles()
    call flux()
    call solve()
    continue_solve = 0
  end do
  continue_timestep = 0
end do

call write_output()
call do_finalisation()
end program
