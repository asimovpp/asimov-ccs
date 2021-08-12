submodule (flowsolve) stepper_3loop
contains
  module subroutine stepper()
    integer :: timestep, solvestep, extrastep

    print *,"Using 3 loop stepper"
    do timestep = 1, 3 
      call do_time_step_stuff()
      do solvestep = 1, 2
        call turbulence()
        ! call particles()
        do extrastep = 1, 2
          call flux()
          call solve()
        end do  
      end do
    end do
  end subroutine stepper
end submodule stepper_3loop

