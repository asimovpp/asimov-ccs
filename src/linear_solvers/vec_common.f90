submodule (vec) vec_common

  use constants, only: cell
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
    class(ccs_mesh), target, intent(in) :: geometry
    type(vector_init_data), intent(inout) :: vector_descriptor

    vector_descriptor%par_env => par_env
    vector_descriptor%mesh => geometry
  end subroutine set_vector_size

  !> @brief Set vector values to be located at either cell-centre or face
  !
  module subroutine set_vector_location(loc, vector_descriptor)
    integer(ccs_int), intent(in) :: loc
    type(vector_init_data), intent (inout) :: vector_descriptor

    vector_descriptor%loc = loc
  end subroutine set_vector_location

  module procedure create_vector_values
    allocate(val_dat%idx(nrows))
    allocate(val_dat%val(nrows))
  end procedure create_vector_values

  module procedure set_vector_values_mode
    val_dat%mode = mode
  end procedure set_vector_values_mode
  
  module procedure set_vector_values_entry

    use constants, only : add_mode, insert_mode
    
    associate(x => val_dat%val(val_dat%current_entry), &
         mode => val_dat%mode)
      if (mode == insert_mode) then
        x = val
      else if (mode == add_mode) then
        x = x + val
      else
        print *, "ERROR: Unrecognised entry mode ", mode
        stop
      end if
    end associate
    
  end procedure set_vector_values_entry
  
end submodule 
