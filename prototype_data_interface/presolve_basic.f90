submodule (presolve) basic_presolve
contains
  
  module subroutine do_initialisation(all_d)
    type(all_data), intent(inout) :: all_d
    print *,"Initialising basic"
    call initialise_turbulence(all_d)
  end subroutine do_initialisation



end submodule basic_presolve
