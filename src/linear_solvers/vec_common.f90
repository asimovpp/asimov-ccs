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

end submodule