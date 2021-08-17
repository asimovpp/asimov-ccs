submodule (solver) solve_cgstab
implicit none
contains
  module subroutine initialise_solve(d)
    type(all_data), intent(inout) :: d
    
    allocate(solver_cgstab_data :: d%solver_d)

    select type (t => d%solver_d)
      type is (solver_cgstab_data)
        print *,"Initialising cgstab solver"
        t%x1 = 77777
    end select
  end subroutine initialise_solve

  module subroutine solve(d)
    class(solver_data), intent(inout) :: d
    type(solver_cgstab_data), pointer :: my_data
    select type (d)
      type is (solver_cgstab_data)
        my_data => d
      class default
        print *, "Unexpected type in cgstab solve"
        stop 1
    end select

    print *,"Solving CG-STAB", my_data%x1
  end subroutine solve
end submodule solve_cgstab
