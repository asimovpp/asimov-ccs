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
  end subroutine initialise_vector

  !> @brief Setter for vector size
  !
  !> param[in/out] vector_descriptor - the vector data object
  !> param[in] mesh                  - the mesh - contains the
  !!                                   information to set the
  !!                                   vector size
  !> param[in] par_env               - the parallel environment 
  !!                                   where the vector resides
  module procedure set_vector_size
    vector_descriptor%par_env => par_env
    vector_descriptor%mesh => geometry
  end procedure

end submodule 
