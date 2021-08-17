submodule (solver) solve_amg
implicit none
contains
  module subroutine initialise_solve(d)
    type(all_data), intent(inout) :: d
    
    allocate(solver_amg_data :: d%solver_d)

    select type (t => d%solver_d)
      type is (solver_amg_data)
        print *,"Initialising amg solver"
        t%x1 = 99999

        select type (fpd => d%form_problem_d)
          type is (form_problem_data_basic)
            t%datar => fpd%datar 
        end select
    end select
  end subroutine initialise_solve

  module subroutine solve(d)
    class(solver_data), intent(inout) :: d
    type(solver_amg_data), pointer :: my_data
    select type (d)
      type is (solver_amg_data)
        my_data => d
      class default
        print *, "Unexpected type in amg solve"
        stop 1
    end select

    print *,"Solving AMG", my_data%x1, my_data%datar(1), my_data%datar(4)
  end subroutine solve
end submodule solve_amg
