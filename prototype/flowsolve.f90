module flowsolve
implicit none

interface
  module subroutine turbulance()
  end subroutine turbulance 
  
  module subroutine particles()
  end subroutine particles
  
  module subroutine flux()
  end subroutine flux
  
  module subroutine solve()
  end subroutine solve
end interface

contains
  subroutine do_time_step_stuff()
    print *,"Doing timestep"
  end subroutine do_time_step_stuff

end module flowsolve
