submodule (flowsolve) stepper_2loop
contains
  module subroutine stepper()
    integer :: timestep, solvestep

    print *,"Using 2 loop stepper"
    do timestep = 1, 3 
      call do_time_step_stuff()
      do solvestep = 1, 2
        call turbulence()
        call particles()
        call flux()
        call solve()
      end do
    end do
  end subroutine stepper
end submodule stepper_2loop

