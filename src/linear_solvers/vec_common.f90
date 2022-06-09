submodule (vec) vec_common
#include "ccs_macros.inc"

  use utils, only : exit_print, str
  use constants, only: cell
  implicit none 

contains

  !>  Constructor for default vector values
  module subroutine initialise_vector(vec_properties)
    type(vector_spec), intent(inout) :: vec_properties    !< the initialised vector values
    vec_properties%par_env => null()
    vec_properties%mesh => null()
    vec_properties%storage_location = cell  ! Default to cell-centre values (so as not to break previous work)
  end subroutine initialise_vector

  !>  Setter for vector size
  module subroutine set_vector_size(par_env, mesh, vec_properties)
    class(parallel_environment), allocatable, target, intent(in) :: par_env   !< the parallel environment where the vector resides
    class(ccs_mesh), target, intent(in) :: mesh                               !< the mesh - contains the information to set the vector size
    type(vector_spec), intent(inout) :: vec_properties                        !< the vector data object

    vec_properties%par_env => par_env
    vec_properties%mesh => mesh
  end subroutine set_vector_size

  !>  Set vector values to be located at either cell-centre or face
  !
  module subroutine set_vector_location(loc, vec_properties)
    integer(ccs_int), intent(in) :: loc
    type(vector_spec), intent (inout) :: vec_properties

    vec_properties%storage_location = loc
  end subroutine set_vector_location

  module procedure create_vector_values
    allocate(val_dat%global_indices(nrows))
    allocate(val_dat%values(nrows))

    val_dat%global_indices(:) = -1_ccs_int
    val_dat%values(:) = 0.0_ccs_real
  end procedure create_vector_values

  module procedure set_vector_values_mode
    val_dat%setter_mode = mode
  end procedure set_vector_values_mode
  
  module procedure set_vector_values_entry

    use constants, only : add_mode, insert_mode
    
    associate(x => val_dat%values(val_dat%current_entry), &
         mode => val_dat%setter_mode)
      if (mode == insert_mode) then
        x = val
      else if (mode == add_mode) then
        x = x + val
      else
        call error_abort("ERROR: Unrecognised entry mode " // str(mode))
      end if
    end associate
    
  end procedure set_vector_values_entry
  
end submodule 
