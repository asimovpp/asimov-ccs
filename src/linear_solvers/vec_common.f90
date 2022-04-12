submodule (vec) vec_common

  use constants, only: cell
  implicit none 

contains

  !> @brief Constructor for default vector values
  !
  !> param[in/out] vec_properties - the initialised vector values
  module subroutine initialise_vector(vec_properties)
    type(vector_spec), intent(inout) :: vec_properties
    vec_properties%par_env => null()
    vec_properties%mesh => null()
    vec_properties%storage_location = cell  ! Default to cell-centre values (so as not to break previous work)
  end subroutine initialise_vector

  !> @brief Setter for vector size
  !
  !> param[in] par_env               - the parallel environment 
  !!                                   where the vector resides
  !> param[in] mesh                  - the mesh - contains the
  !!                                   information to set the
  !!                                   vector size
  !> param[in/out] vec_properties - the vector data object
  module subroutine set_vector_size(par_env, geometry, vec_properties)
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    class(ccs_mesh), target, intent(in) :: geometry
    type(vector_spec), intent(inout) :: vec_properties

    vec_properties%par_env => par_env
    vec_properties%mesh => geometry
  end subroutine set_vector_size

  !> @brief Set vector values to be located at either cell-centre or face
  !
  module subroutine set_vector_location(loc, vec_properties)
    integer(ccs_int), intent(in) :: loc
    type(vector_spec), intent (inout) :: vec_properties

    vec_properties%storage_location = loc
  end subroutine set_vector_location

end submodule 
