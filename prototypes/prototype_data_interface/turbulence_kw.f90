submodule (flowsolve) turbulence_kw
contains
  module subroutine turbulence(input)
    class(turbulence_input), intent(inout) :: input
    select type (input)
      type is (turb_2)
        print *,"Calculating turbulence kw", input%c2, input%c3
      class default
        print *, "Unexpected type in turbulence kw"
        stop 1
    end select
  end subroutine turbulence
end submodule turbulence_kw

submodule (presolve) initialise_turbulence_kw
contains
  module subroutine initialise_turbulence(all_d)
    type(all_data), intent(inout) :: all_d
    allocate(turb_2 :: all_d%turb_data)
    select type (t => all_d%turb_data)
      type is (turb_2)
        t%c2 = 42
        t%c3 = 99
    end select
    print *,"Initialising turbulence kw"
  end subroutine initialise_turbulence
end submodule initialise_turbulence_kw
