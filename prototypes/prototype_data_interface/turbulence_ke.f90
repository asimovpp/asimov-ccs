submodule (flowsolve) turbulence_ke
contains
  module subroutine turbulence(input)
    class(turbulence_input), intent(inout) :: input
    type(turb_1), pointer :: my_data
    select type (input)
      type is (turb_1)
        my_data => input
      class default
        print *, "Unexpected type in turbulence ke"
        stop 1
    end select

    print *,"Calculating turbulence ke", my_data%c1

  end subroutine turbulence
end submodule turbulence_ke 


submodule (presolve) initialise_turbulence_ke
contains
  module subroutine initialise_turbulence(all_d)
    type(all_data), intent(inout) :: all_d

    allocate(turb_1 :: all_d%turb_data)

    select type (t => all_d%turb_data)
      type is (turb_1)
        t%c1 = 33
    end select
    print *,"Initialising turbulence ke"
  end subroutine initialise_turbulence
end submodule initialise_turbulence_ke
