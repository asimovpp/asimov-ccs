submodule (mat) mat_common
#include "ccs_macros.inc"

  use utils, only : exit_print, str

  implicit none

contains

  !>  Constructor for default matrix values
  module subroutine initialise_matrix(mat_properties)
    type(matrix_spec), intent(inout) :: mat_properties  !< the initialised matrix values
    mat_properties%par_env => null()
    mat_properties%par_env => null()
  end subroutine initialise_matrix

  !>  Setter for global matrix size
  module subroutine set_matrix_size(par_env, mesh, mat_properties)
    class(parallel_environment), allocatable, target, intent(in) :: par_env   !< the parallel environment where the matrix resides
    class(ccs_mesh), target, intent(in) :: mesh                               !< the mesh object
    type(matrix_spec), intent(inout) :: mat_properties                        !< the matrix data object

    mat_properties%mesh => mesh
    mat_properties%par_env => par_env
  end subroutine set_matrix_size

  !>  Setter for matrix number of non-zeros
  module subroutine set_nnz(nnz, mat_properties)
    integer(ccs_int), intent(in) :: nnz                 !< the number of non-zeros
    type(matrix_spec), intent(inout) :: mat_properties  !< the matrix data object

    mat_properties%nnz = nnz
  end subroutine  set_nnz

  module procedure create_matrix_values
    allocate(val_dat%global_row_indices(nrows))
    allocate(val_dat%global_col_indices(nrows))
    allocate(val_dat%values(nrows))
  end procedure create_matrix_values

  module procedure set_matrix_values_mode
    val_dat%setter_mode = mode
  end procedure set_matrix_values_mode
  
  module subroutine set_matrix_values_entry(val, val_dat)

    use constants, only : add_mode, insert_mode

    real(ccs_real), intent(in) :: val
    type(matrix_values), intent(inout) :: val_dat

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
    
  end subroutine set_matrix_values_entry
  
end submodule 
