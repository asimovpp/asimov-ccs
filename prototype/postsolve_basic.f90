submodule (postsolve) basic_postsolve
contains
  module subroutine write_output()
    print *,"Writing output basic"
  end subroutine write_output
  
  module subroutine do_finalisation()
    print *,"Finalising basic"
  end subroutine do_finalisation
end submodule basic_postsolve
