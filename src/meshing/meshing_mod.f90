!> @brief Module file meshing.mod
!>
!> @details Module defining meshing interface for ASiMoV-CCS

module meshing

  use constants, only : ndim
  use kinds, only : ccs_int, ccs_real
  use types, only : mesh, face_locator, cell_locator, neighbour_locator
  
  implicit none

  private
  public :: set_face_location
  public :: set_cell_location
  public :: set_neighbour_location
  public :: set_face_index
  public :: get_face_normal
  public :: get_face_area
  public :: get_centre
  public :: get_volume
  public :: get_global_index
  public :: get_local_index
  public :: get_boundary_status
  public :: get_local_status
  public :: count_neighbours
  public :: get_distance
  
  interface get_centre
    module procedure get_cell_centre
    module procedure get_neighbour_centre
    module procedure get_face_centre
  end interface get_centre

  interface get_global_index
    module procedure get_cell_global_index
    module procedure get_neighbour_global_index
  end interface get_global_index

  interface get_local_index
    module procedure get_cell_local_index
    module procedure get_neighbour_local_index
    module procedure get_face_local_index
  end interface get_local_index
  
  interface count_neighbours
    module procedure cell_count_neighbours
  end interface count_neighbours

  interface get_boundary_status
    module procedure get_neighbour_boundary_status
    module procedure get_face_boundary_status
  end interface get_boundary_status

  interface get_volume
    module procedure get_cell_volume
    module procedure get_neighbour_volume
  end interface get_volume
  
  interface get_distance
    module procedure get_neighbour_distance
    module procedure get_face_distance
  end interface get_distance
  
  interface

    !> @brief Constructs a cell locator object.
    !
    !> @description Creates the association between a mesh and cell index, storing it in the
    !!              returned cell locator object.
    !
    !> @param[in]  mesh         geometry      - the mesh object being referred to.
    !> @param[in]  ccs_int     cell_idx      - the cell index. 
    !> @param[out] cell_locator cell_location - the cell locator object linking a cell index with
    !!                                          the mesh.
  module subroutine set_cell_location(geometry, cell_idx, cell_location)
      type(mesh), target, intent(in) :: geometry
      integer(ccs_int), intent(in) :: cell_idx
      type(cell_locator), intent(out) :: cell_location
    end subroutine set_cell_location
  
    !> @brief Constructs a face locator object.
    !
    !> @description Creates the association between a face relative to a cell, i.e. to access the
    !!              nth face of cell i.
    !
    !> @param[in]  mesh         geometry      - the mesh object being referred to.
    !> @param[in]  ccs_int     cell_idx      - the index of the cell whose face is being accessed.
    !> @param[in]  ccs_int     cell_face_ctr - the cell-local index of the face.
    !> @param[out] face_locator face_location - the face locator object linking a cell-relative
    !!                                          index with the mesh.
    module subroutine set_face_location(geometry, cell_idx, cell_face_ctr, face_location)
      type(mesh), target, intent(in) :: geometry
      integer(ccs_int), intent(in) :: cell_idx
      integer(ccs_int), intent(in) :: cell_face_ctr
      type(face_locator), intent(out) :: face_location
    end subroutine set_face_location

    !> @brief Constructs a neighbour locator object.
    !
    !> @description Creates the association between a neighbour cell F relative to cell P, i.e. to
    !!              access the nth neighbour of cell i.
    !
    !> @param[in]  cell_locator      cell_location      - the cell locator object of the cell whose
    !!                                                    neighbour is being accessed.
    !> @param[in]  ccs_int          cell_neighbour_ctr - the cell-local index of the neighbour.
    !> @param[out] neighbour_locator neighbour_location - the neighbour locator object linking a
    !!                                                    cell-relative index with the mesh.
    module subroutine set_neighbour_location(cell_location, cell_neighbour_ctr, neighbour_location)
      type(cell_locator), intent(in) :: cell_location
      integer(ccs_int), intent(in) :: cell_neighbour_ctr
      type(neighbour_locator), intent(out) :: neighbour_location
    end subroutine set_neighbour_location

    !> @brief Set face index
    module subroutine set_face_index(cell_idx, cell_face_ctr, face_idx, geometry)
      integer(ccs_int), intent(in) :: cell_idx
      integer(ccs_int), intent(in) :: cell_face_ctr
      integer(ccs_int), intent(in) :: face_idx
      type(mesh), target, intent(inout) :: geometry
    end subroutine set_face_index

    !> @brief Returns the normal vector of a face
    !
    !> @param[in]  face_locator    face_location - the face locator object.
    !> @param[out] real(ccs_real) normal(ndim)  - an ndimensional array representing the face normal
    !!                                             vector.
    module subroutine get_face_normal(face_location, normal)
      type(face_locator), intent(in) :: face_location
      real(ccs_real), dimension(ndim), intent(out) :: normal
    end subroutine get_face_normal

    !> @brief Returns the area of a face
    !
    !> @param[in]  face_locator    face_location - the face locator object.
    !> @param[out] real(ccs_real) area          - the face area.
    module subroutine get_face_area(face_location, area)
      type(face_locator), intent(in) :: face_location
      real(ccs_real), intent(out) :: area
    end subroutine get_face_area

    !> @brief Returns the centre of a cell
    !
    !> @param[in]  cell_locator     cell_location - the cell locator object.
    !> @param[out] real(ccs_real)  x(ndim)       - an ndimensional array representing the cell centre.
    module subroutine get_cell_centre(cell_location, x)
      type(cell_locator), intent(in) :: cell_location
      real(ccs_real), dimension(ndim), intent(out) :: x
    end subroutine get_cell_centre

    !> @brief Returns the centre of a neighbour cell
    !
    !> @param[in]  neighbour_locator neighbour_location - the neighbour locator object.
    !> @param[out] real(ccs_real)   x(ndim)            - an ndimensional array representing the
    !!                                                    neighbour cell centre.
    module subroutine get_neighbour_centre(neighbour_location, x)
      type(neighbour_locator), intent(in) :: neighbour_location
      real(ccs_real), dimension(ndim), intent(out) :: x
    end subroutine get_neighbour_centre

    !> @brief Returns the centre of a face
    !
    !> @param[in]  face_locator     face_location - the face locator object.
    !> @param[out] real(ccs_real)  x(ndim)       - an ndimensional array representing the face centre.
    module subroutine get_face_centre(face_location, x)
      type(face_locator), intent(in) :: face_location
      real(ccs_real), dimension(ndim), intent(out) :: x
    end subroutine get_face_centre

    !> @brief Returns the volume of a cell
    !
    !> @param[in] cell_locator     cell_location - the cell locator object.
    !> @param[out] real(ccs_real) V             - the cell volume.
    module subroutine get_cell_volume(cell_location, V)
      type(cell_locator), intent(in) :: cell_location
      real(ccs_real), intent(out) :: V
    end subroutine get_cell_volume

    !> @brief Returns the volume of a neighbour cell
    !
    !> @param[in] neighbour_locator neighbour_location - the neighbour locator object.
    !> @param[out] real(ccs_real)  V                  - the neighbour cell volume.
    module subroutine get_neighbour_volume(neighbour_location, V)
      type(neighbour_locator), intent(in) :: neighbour_location
      real(ccs_real), intent(out) :: V
    end subroutine get_neighbour_volume

    !> @brief Returns the global index of a cell
    !
    !> @param[in]  cell_locator      cell_location - the cell locator object.
    !> @param[out] integer(ccs_int) idxg          - the global index of the cell.
    module subroutine get_cell_global_index(cell_location, idxg)
      type(cell_locator), intent(in) :: cell_location
      integer(ccs_int), intent(out) :: idxg
    end subroutine get_cell_global_index

    !> @brief Returns the global index of a neighbouring cell
    !
    !> @param[in]  neighbour_locator neighbour_location - the neighbour locator object.
    !> @param[out] integer(ccs_int) nbidxg             - the global index of the neighbour cell.
    module subroutine get_neighbour_global_index(neighbour_location, nbidxg)
      type(neighbour_locator), intent(in) :: neighbour_location
      integer(ccs_int), intent(out) :: nbidxg
    end subroutine get_neighbour_global_index

    !> @brief Returns the neighbour count of a cell (including boundary neighbours)
    !
    !> @param[in]  cell_locator      cell_location - the cell locator object.
    !> @param[out] integer(ccs_int) nnb           - the neighbour count of the cell.
    module subroutine cell_count_neighbours(cell_location, nnb)
      type(cell_locator), intent(in) :: cell_location
      integer(ccs_int), intent(out) :: nnb
    end subroutine cell_count_neighbours

    !> @brief Returns the boundary status of a neighbouring cell
    !
    !> @param[in]  neighbour_locator neighbour_location - the neighbour locator object.
    !> @param[out] logical           is_boundary        - the boundary status of the neighbour.
    module subroutine get_neighbour_boundary_status(neighbour_location, is_boundary)
      type(neighbour_locator), intent(in) :: neighbour_location
      logical, intent(out) :: is_boundary
    end subroutine get_neighbour_boundary_status

    !> @brief Returns the boundary status of a face
    !
    !> @param[in]  face_locator face_location - the face locator object.
    !> @param[out] logical      is_boundary   - the boundary status of the neighbour.
    module subroutine get_face_boundary_status(face_location, is_boundary)
      type(face_locator), intent(in) :: face_location
      logical, intent(out) :: is_boundary
    end subroutine get_face_boundary_status

    !> @brief Returns the local distribution status of a neighbouring cell
    !
    !> @description Given a distributed mesh, a processor needs both the cells within its partition
    !!              and cells from the surrounding halo - this subroutine get_indicates whether a
    !!              cell's neighbour is within the local partition or the halo.
    !
    !> @param[in]  neighbour_locator neighbour_location - the neighbour locator object.
    !> @param[out] logical           is_local           - the local status of the neighbour.
    module subroutine get_local_status(neighbour_location, is_local)
      type(neighbour_locator), intent(in) :: neighbour_location
      logical, intent(out) :: is_local
    end subroutine get_local_status

    !> @brief Returns the local index of a cell
    !
    !> @description Generally the local index of a cell is should be the same as its location within
    !!              the local cell vector - this particular subroutine get_is therefore expected of
    !!              limited use and is mostly present for uniformity.
    !
    !> @param[in]  cell_locator      cell_location - the cell locator object.
    !> @param[out] integer(ccs_int) idx           - the local index of the cell.
    module subroutine get_cell_local_index(cell_location, idx)
      type(cell_locator), intent(in) :: cell_location
      integer(ccs_int), intent(out) :: idx
    end subroutine get_cell_local_index

    !> @brief Returns the local index of a neighbouring cell
    !
    !> @param[in]  neighbour_locator neighbour_location - the neighbour locator object.
    !> @param[out] integer(ccs_int) nbidx              - the local index of the neighbour cell.
    module subroutine get_neighbour_local_index(neighbour_location, nbidx)
      type(neighbour_locator), intent(in) :: neighbour_location
      integer(ccs_int), intent(out) :: nbidx
    end subroutine get_neighbour_local_index

    !> @brief Returns the local index of a face
    !
    !> @param[in]  face_locator      face_location - the face locator object.
    !> @param[out] integer(ccs_int) idx           - the local index of the face.
    module subroutine get_face_local_index(face_location, idx)
      type(face_locator), intent(in) :: face_location
      integer(ccs_int), intent(out) :: idx
    end subroutine get_face_local_index

    !> @brief Returns the distance between two cell centres
    !
    !> @param[in]  cell_locator      loc_p  - The cell distance is measured from.
    !> @param[in]  neighbour_locator loc_nb - The cell distance is measured to.
    !> @param[out] ccs_real         dx     - ndim-array of the distance
    module subroutine get_neighbour_distance(loc_p, loc_nb, dx)
      type(cell_locator), intent(in) :: loc_p
      type(neighbour_locator), intent(in) :: loc_nb
      real(ccs_real), dimension(ndim), intent(out) :: dx
    end subroutine get_neighbour_distance

    !> @brief Returns the distance from cell to face centres
    !
    !> @param[in]  cell_locator loc_p - The cell distance is measured from.
    !> @param[in]  face_locator loc_f - The face distance is measured to.
    !> @param[out] ccs_real    dx    - ndim-array of the distance
    module subroutine get_face_distance(loc_p, loc_f, dx)
      type(cell_locator), intent(in) :: loc_p
      type(face_locator), intent(in) :: loc_f
      real(ccs_real), dimension(ndim), intent(out) :: dx
    end subroutine get_face_distance
  end interface

end module meshing
  
