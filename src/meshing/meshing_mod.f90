!> @brief Module file meshing.mod
!>
!> @details Module defining meshing interface for ASiMoV-CCS

module meshing

  use constants, only : ndim
  use kinds, only : accs_int, accs_real
  use types, only : mesh, face_locator, cell_locator, neighbour_locator
  
  implicit none

  private
  public :: set_face_location
  public :: set_cell_location
  public :: set_neighbour_location
  public :: face_normal
  public :: face_area
  public :: centre
  public :: volume
  public :: global_index
  public :: local_index
  public :: count_neighbours
  public :: boundary_status
  public :: local_status

  interface centre
    module procedure cell_centre
    module procedure face_centre
  end interface centre

  interface global_index
    module procedure cell_global_index
    module procedure neighbour_global_index
  end interface global_index

  interface local_index
    module procedure cell_local_index
    module procedure neighbour_local_index
  end interface local_index
  
  interface count_neighbours
    module procedure cell_count_neighbours
  end interface count_neighbours
  
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

    !> @brief Returns the normal vector of a face
    !
    !> @param[in]  face_locator    face_location - the face locator object.
    !> @param[out] real(accs_real) normal(ndim)  - an ndimensional array representing the face normal
    !!                                             vector.
    module subroutine face_normal(face_location, normal)
      type(face_locator), intent(in) :: face_location
      real(accs_real), dimension(ndim), intent(out) :: normal
    end subroutine face_normal

    !> @brief Returns the area of a face
    !
    !> @param[in]  face_locator    face_location - the face locator object.
    !> @param[out] real(accs_real) area          - the face area.
    module subroutine face_area(face_location, area)
      type(face_locator), intent(in) :: face_location
      real(accs_real), intent(out) :: area
    end subroutine face_area

    !> @brief Returns the centre of a cell
    !
    !> @param[in]  cell_locator     cell_location - the cell locator object.
    !> @param[out] real(accs_real)  x(ndim)       - an ndimensional array representing the cell centre.
    module subroutine cell_centre(cell_location, x)
      type(cell_locator), intent(in) :: cell_location
      real(accs_real), dimension(ndim), intent(out) :: x
    end subroutine cell_centre

    !> @brief Returns the centre of a face
    !
    !> @param[in]  face_locator     face_location - the face locator object.
    !> @param[out] real(accs_real)  x(ndim)       - an ndimensional array representing the face centre.
    module subroutine face_centre(face_location, x)
      type(face_locator), intent(in) :: face_location
      real(accs_real), dimension(ndim), intent(out) :: x
    end subroutine face_centre

    !> @brief Returns the volume of a cell
    !
    !> @param[in] cell_locator     cell_location - the cell locator object.
    !> @param[out] real(accs_real) V             - the cell volume.
    module subroutine volume(cell_location, V)
      type(cell_locator), intent(in) :: cell_location
      real(accs_real), intent(out) :: V
    end subroutine volume

    !> @brief Returns the global index of a cell
    !
    !> @param[in]  cell_locator      cell_location - the cell locator object.
    !> @param[out] integer(accs_int) idxg          - the global index of the cell.
    module subroutine cell_global_index(cell_location, idxg)
      type(cell_locator), intent(in) :: cell_location
      integer(accs_int), intent(out) :: idxg
    end subroutine cell_global_index

    !> @brief Returns the global index of a neighbouring cell
    !
    !> @param[in]  neighbour_locator neighbour_location - the neighbour locator object.
    !> @param[out] integer(accs_int) nbidxg             - the global index of the neighbour cell.
    module subroutine neighbour_global_index(neighbour_location, nbidxg)
      type(neighbour_locator), intent(in) :: neighbour_location
      integer(accs_int), intent(out) :: nbidxg
    end subroutine neighbour_global_index

    !> @brief Returns the neighbour count of a cell (including boundary neighbours)
    !
    !> @param[in]  cell_locator      cell_location - the cell locator object.
    !> @param[out] integer(accs_int) nnb           - the neighbour count of the cell.
    module subroutine cell_count_neighbours(cell_location, nnb)
      type(cell_locator), intent(in) :: cell_location
      integer(accs_int), intent(out) :: nnb
    end subroutine cell_count_neighbours

    !> @brief Returns the boundary status of a neighbouring cell
    !
    !> @param[in]  neighbour_locator neighbour_location - the neighbour locator object.
    !> @param[out] logical           is_boundary        - the boundary status of the neighbour.
    module subroutine boundary_status(neighbour_location, is_boundary)
      type(neighbour_locator), intent(in) :: neighbour_location
      logical, intent(out) :: is_boundary
    end subroutine boundary_status

    !> @brief Returns the local distribution status of a neighbouring cell
    !
    !> @description Given a distributed mesh, a processor needs both the cells within its partition
    !!              and cells from the surrounding halo - this subroutine indicates whether a
    !!              cell's neighbour is within the local partition or the halo.
    !
    !> @param[in]  neighbour_locator neighbour_location - the neighbour locator object.
    !> @param[out] logical           is_local           - the local status of the neighbour.
    module subroutine local_status(neighbour_location, is_local)
      type(neighbour_locator), intent(in) :: neighbour_location
      logical, intent(out) :: is_local
    end subroutine local_status

    !> @brief Returns the local index of a cell
    !
    !> @description Generally the local index of a cell is should be the same as its location within
    !!              the local cell vector - this particular subroutine is therefore expected of
    !!              limited use and is mostly present for uniformity.
    !
    !> @param[in]  cell_locator      cell_location - the cell locator object.
    !> @param[out] integer(accs_int) idx           - the local index of the cell.
    module subroutine cell_local_index(cell_location, idx)
      type(cell_locator), intent(in) :: cell_location
      integer(accs_int), intent(out) :: idx
    end subroutine cell_local_index

    !> @brief Returns the local index of a neighbouring cell
    !
    !> @param[in]  neighbour_locator neighbour_location - the neighbour locator object.
    !> @param[out] integer(accs_int) nbidx              - the local index of the neighbour cell.
    module subroutine neighbour_local_index(neighbour_location, nbidx)
      type(neighbour_locator), intent(in) :: neighbour_location
      integer(accs_int), intent(out) :: nbidx
    end subroutine neighbour_local_index
  end interface

end module meshing
  
