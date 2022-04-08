submodule (vec) vec_common

  implicit none 

contains

  !> @brief Constructor for default vector values
  !
  !> param[in/out] vector_descriptor - the initialised vector values
  module subroutine initialise_vector(vector_descriptor)
    type(vector_init_data), intent(inout) :: vector_descriptor
    vector_descriptor%par_env => null()
    vector_descriptor%mesh => null()
    vector_descriptor%loc = cell  ! Default to cell-centre values (so as not to break previous work)
  end subroutine initialise_vector

  !> @brief Setter for vector size
  !
  !> param[in] par_env               - the parallel environment 
  !!                                   where the vector resides
  !> param[in] mesh                  - the mesh - contains the
  !!                                   information to set the
  !!                                   vector size
  !> param[in/out] vector_descriptor - the vector data object
  module subroutine set_vector_size(par_env, geometry, vector_descriptor)
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    class(mesh), target, intent(in) :: geometry
    type(vector_init_data), intent(inout) :: vector_descriptor

    vector_descriptor%par_env => par_env
    vector_descriptor%mesh => geometry
  end subroutine set_vector_size

  !> @brief Set vector values to be located at either cell-centre or face
  !
  module subroutine set_vector_location(loc, vector_descriptor)
    integer(accs_int), intent(in) :: loc
    type(vector_init_data), intent (inout) :: vector_descriptor

    vector_descriptor%loc = loc
  end subroutine set_vector_location

end submodule 
