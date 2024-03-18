submodule(reordering) reordering_none
#include "ccs_macros.inc"

  use utils, only: debug_print

  implicit none

contains

  !v Creates a no-action reordering list
  module subroutine get_reordering(mesh, new_indices)

    use meshing, only: get_local_num_cells

    type(ccs_mesh), intent(in) :: mesh                                      !< the mesh to be reordered
    integer(ccs_int), dimension(:), allocatable, intent(out) :: new_indices !< new indices in "to(from)" format

    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: i

    call dprint("Reordering with no-reordering.")
    call get_local_num_cells(mesh, local_num_cells)

    allocate (new_indices(local_num_cells))
    new_indices = (/(i, i=1,local_num_cells)/)

  end subroutine get_reordering

end submodule
