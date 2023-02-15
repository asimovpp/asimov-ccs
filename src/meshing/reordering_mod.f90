module reordering

  use kinds, only: ccs_int
  use types, only: ccs_mesh

  implicit none

  private
  public :: reorder_cells
  public :: bandwidth

interface
  
  module subroutine reorder_cells(mesh)
    type(ccs_mesh), intent(inout) :: mesh
  end subroutine
  
  module subroutine apply_reordering(new_indices, mesh)
    integer(ccs_int), dimension(:), intent(in) :: new_indices
    type(ccs_mesh), intent(inout) :: mesh
  end subroutine
  
  module subroutine bandwidth(mesh)
    type(ccs_mesh), intent(in) :: mesh
  end subroutine
  
  module subroutine get_reordering(mesh, new_indices)
    type(ccs_mesh), intent(in) :: mesh
    integer(ccs_int), dimension(:), allocatable, intent(out) :: new_indices
  end subroutine
  

end interface

end module reordering
