submodule (turbulence) turbulence_kw
implicit none
contains

  module subroutine initialise_turbulence(d)
    type(all_data), intent(inout) :: d
    allocate(turbulence_data_2 :: d%turbulence_d)
    select type (t => d%turbulence_d)
      type is (turbulence_data_2)
        t%c2 = 42
        t%c3 = 99
    end select
    print *,"Initialising turbulence kw"
  end subroutine initialise_turbulence

  module subroutine do_turbulence(d)
    class(turbulence_data), intent(inout) :: d 
    select type (d)
      type is (turbulence_data_2)
        print *,"Calculating turbulence kw", d%c2, d%c3
      class default
        print *, "Unexpected type in turbulence kw"
        stop 1
    end select
  end subroutine do_turbulence

end submodule turbulence_kw
