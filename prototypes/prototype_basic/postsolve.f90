module postsolve
implicit none

interface
  module subroutine write_output()
  end subroutine write_output
  
  module subroutine do_finalisation()
  end subroutine do_finalisation
end interface

end module postsolve
