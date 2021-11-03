submodule (vec) vec_common

  implicit none 

contains

  !> @brief Constructor for default vector values
  !
  !> param[in/out] vector_descriptor - the initialised vector values
  module subroutine initialise_vector(vector_descriptor)
    type(vector_init_data), intent(inout) :: vector_descriptor
    vector_descriptor%nglob = -1
    vector_descriptor%nloc = -1
    vector_descriptor%par_env => null()
  end subroutine initialise_vector

  !> @brief Setter for global vector size
  !
  !> param[in/out] vector_descriptor - the vector data object
  !> param[in] size                  - the global vector size
  !> param[in] par_env               - the parallel environment 
  !!                                   where the vector resides
  module subroutine set_global_vector_size(vector_descriptor, size, par_env)
    type(vector_init_data), intent(inout) :: vector_descriptor
    integer(accs_int), intent(in) :: size
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    vector_descriptor%nglob = size
    vector_descriptor%nloc = -1
    vector_descriptor%par_env => par_env
  end subroutine

  !> @brief Setter for local vector size
  !
  !> param[in/out] vvector_descriptor - the vector data object
  !> param[in] size                   - the local vector size
  !> param[in] par_env                - the parallel environment 
  !!                                    where the vector resides
  module subroutine set_local_vector_size(vector_descriptor, size, par_env)
    type(vector_init_data), intent(inout) :: vector_descriptor
    integer(accs_int), intent(in) :: size
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    vector_descriptor%nloc = size
    vector_descriptor%nglob = -1
    vector_descriptor%par_env => par_env
  end subroutine

end submodule 