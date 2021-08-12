submodule (turbulence) turbulence_ke
implicit none
contains
  
  module subroutine initialise_turbulence(d)
    type(all_data), intent(inout) :: d

    allocate(turbulence_data_1 :: d%turbulence_d)

    select type (t => d%turbulence_d)
      type is (turbulence_data_1)
        t%c1 = 33
    end select
    print *,"Initialising turbulence ke"
  end subroutine initialise_turbulence
  
  module subroutine do_turbulence(d)
    class(turbulence_data), intent(inout) :: d
    type(turbulence_data_1), pointer :: my_data
    select type (d)
      type is (turbulence_data_1)
        my_data => d
      class default
        print *, "Unexpected type in turbulence ke"
        stop 1
    end select

    print *,"Calculating turbulence ke", my_data%c1

  end subroutine do_turbulence

end submodule turbulence_ke
