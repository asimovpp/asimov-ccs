submodule (problem_former) form_problem_1
implicit none

contains

  module subroutine initialise_form_problem(d)
    type(all_data), intent(inout) :: d

    allocate(form_problem_data_basic :: d%form_problem_d)

    select type (t => d%form_problem_d)
      type is (form_problem_data_basic)
        t%x1 = 173
        allocate( t%datar(4) )
        t%datar = (/2, 5, 9, 11/)
    end select
    print *,"Initialising basic problem former"
    
  end subroutine initialise_form_problem


  module subroutine form_problem(d)
    class(form_problem_data), intent(inout) :: d
    type(form_problem_data_basic), pointer :: my_data
    select type (d)
      type is (form_problem_data_basic)
        my_data => d
      class default
        print *, "Unexpected type in basic problem former"
        stop 1
    end select

    print *,"Forming problem basic", my_data%x1, my_data%datar(1), my_data%datar(4)
  end subroutine form_problem

end submodule form_problem_1
