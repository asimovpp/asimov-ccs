submodule (vec) vec_common

  implicit none 

contains

  !> @brief Constructor for default vector values
  !
  !> param[in/out] vec     - the initialised vector values
  module subroutine initialise_vector(vec)
    type(vector_init_data), intent(inout) :: vec
    vec%nglob = -1
    vec%nloc = -1
    vec%par_env => null()
  end subroutine initialise_vector

  module subroutine set_global_vector_size(vec, size, par_env)
    type(vector_init_data), intent(inout) :: vec
    integer(accs_int), intent(in) :: size
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    vec%nglob = size
    vec%par_env => par_env
  end subroutine

  module subroutine set_local_vector_size(vec, size, par_env)
    type(vector_init_data), intent(inout) :: vec
    integer(accs_int), intent(in) :: size
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    vec%nloc = size
    vec%par_env => par_env
  end subroutine

end submodule 