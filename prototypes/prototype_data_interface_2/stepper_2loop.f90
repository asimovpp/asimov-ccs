submodule (stepper) stepper_2loop
use turbulence
use solver
implicit none

contains
  module subroutine initialise_stepper(d)
    type(all_data), intent(inout) :: d
    call initialise_turbulence(d)
    call initialise_solve(d)
  end subroutine initialise_stepper

  module subroutine do_steps(d)
    integer :: timestep, solvestep
    type(all_data), intent(inout) :: d

    print *,"Using 2 loop stepper"
    do timestep = 1, 3 
      do solvestep = 1, 2
        call do_turbulence(d%turbulence_d)
        call solve(d%solver_d)
      end do
    end do
  end subroutine do_steps
end submodule stepper_2loop

