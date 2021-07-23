submodule (flowsolve) stepper_2loop
contains
  module subroutine stepper(all_d)
    integer :: timestep, solvestep
    type(all_data), intent(inout) :: all_d

    print *,"Using 2 loop stepper"
    do timestep = 1, 3 
      do solvestep = 1, 2
        call turbulence(all_d%turb_data)
        call solve()
      end do
    end do
  end subroutine stepper
end submodule stepper_2loop

