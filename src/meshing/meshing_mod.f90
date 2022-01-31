!> @brief Module file meshing.mod
!>
!> @details Module defining meshing interface for ASiMoV-CCS

module meshing

  use kinds, only : accs_int
  use types, only : mesh, face_locator, cell_locator, neighbour_locator
  
  implicit none

  private
  public :: set_face_location
  public :: set_cell_location
  public :: set_neighbour_location

  interface

    !> @brief Constructs a cell locator object.
    !
    !> @description Creates the association between a mesh and cell index, storing it in the
    !!              returned cell locator object.
    !
    !> @param[out] cell_locator cell_location - the cell locator object linking a cell index with
    !!                                          the mesh.
    !> @param[in]  mesh         geometry      - the mesh object being referred to.
    !> @param[in]  accs_int     cell_idx      - the cell index. 
    module subroutine set_cell_location(cell_location, geometry, cell_idx)
      type(cell_locator), intent(out) :: cell_location
      type(mesh), target, intent(in) :: geometry
      integer(accs_int), intent(in) :: cell_idx
    end subroutine set_cell_location

    !> @brief Constructs a face locator object.
    !
    !> @description Creates the association between a face relative to a cell, i.e. to access the
    !!              nth face of cell i.
    !
    !> @param[out] face_locator face_location - the face locator object linking a cell-relative
    !!                                          index with the mesh.
    !> @param[in]  mesh         geometry      - the mesh object being referred to.
    !> @param[in]  accs_int     cell_idx      - the index of the cell whose face is being accessed.
    !> @param[in]  accs_int     cell_face_ctr - the cell-local index of the face.
    module subroutine set_face_location(face_location, geometry, cell_idx, cell_face_ctr)
      type(face_locator), intent(out) :: face_location
      type(mesh), target, intent(in) :: geometry
      integer(accs_int), intent(in) :: cell_idx
      integer(accs_int), intent(in) :: cell_face_ctr
    end subroutine set_face_location

    !> @brief Constructs a neighbour locator object.
    !
    !> @description Creates the association between a neighbour cell F relative to cell P, i.e. to
    !!              access the nth neighbour of cell i.
    !
    !> @param[out] neighbour_locator neighbour_location - the neighbour locator object linking a
    !!                                                    cell-relative index with the mesh.
    !> @param[in]  cell_locator      cell_location      - the cell locator object of the cell whose
    !!                                                    neighbour is being accessed.
    !> @param[in]  accs_int          cell_neighbour_ctr - the cell-local index of the neighbour.
    module subroutine set_neighbour_location(neighbour_location, cell_location, cell_neighbour_ctr)
      type(neighbour_locator), intent(out) :: neighbour_location
      type(cell_locator), intent(in) :: cell_location
      integer(accs_int), intent(in) :: cell_neighbour_ctr
    end subroutine set_neighbour_location
  end interface

contains
  
  module subroutine set_face_location(face_location, geometry, cell_idx, cell_face_ctr)
    type(face_locator), intent(out) :: face_location
    type(mesh), target, intent(in) :: geometry
    integer(accs_int), intent(in) :: cell_idx
    integer(accs_int), intent(in) :: cell_face_ctr

    face_location%mesh => geometry
    face_location%cell_idx = cell_idx
    face_location%cell_face_ctr = cell_face_ctr
  end subroutine set_face_location

  module subroutine set_cell_location(cell_location, geometry, cell_idx)
    type(cell_locator), intent(out) :: cell_location
    type(mesh), target, intent(in) :: geometry
    integer(accs_int), intent(in) :: cell_idx

    ! XXX: Potentially expensive...
    if (cell_idx > size(geometry%idx_global)) then
      print *, "ERROR: trying to access cell I don't own!", cell_idx, geometry%nlocal
      stop
    else
      cell_location%mesh => geometry
      cell_location%cell_idx = cell_idx
    end if

  end subroutine set_cell_location

  module subroutine set_neighbour_location(neighbour_location, cell_location, cell_neighbour_ctr)
    type(neighbour_locator), intent(out) :: neighbour_location
    type(cell_locator), intent(in) :: cell_location
    integer(accs_int), intent(in) :: cell_neighbour_ctr

    ! integer(accs_int) :: nnb
    
    neighbour_location%mesh => cell_location%mesh
    neighbour_location%cell_idx = cell_location%cell_idx

    !! XXX: Safe, but would create a circular dependency...
    !! ! XXX: Potentially expensive...
    !! call count_neighbours(cell_location, nnb)
    !! if (cell_neighbour_ctr > nnb) then
    !!   print *, "ERROR: cell has fewer neighbours than neighbour count requested!"
    !!   stop
    !! else if (cell_neighbour_ctr < 1) then
    !!   print *, "ERROR: cell neighbour counter must be >= 1!"
    !! else
    !!   neighbour_location%cell_neighbour_ctr = cell_neighbour_ctr
    !! end if
    
    neighbour_location%cell_neighbour_ctr = cell_neighbour_ctr

    associate(mymesh => neighbour_location%mesh, &
         i => neighbour_location%cell_idx, &
         j => neighbour_location%cell_neighbour_ctr)
      if (mymesh%nbidx(j, i) == i) then
        print *, "ERROR: trying to set self as neighbour! Cell: ", i, j
      end if
    end associate
    
  end subroutine set_neighbour_location
end module meshing
  
