!>  Module file meshing.mod
!>
!>  Module defining meshing interface for ASiMoV-CCS

module meshing

  use constants, only: ndim
  use kinds, only: ccs_int, ccs_real
  use types, only: ccs_mesh, face_locator, cell_locator, neighbour_locator

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

    !>  Constructs a cell locator object.
    !
    !> @description Creates the association between a mesh and cell index, storing it in the
    !!              returned cell locator object.
    !
    !> @param[in]  mesh         mesh      - the mesh object being referred to.
    !> @param[in]  ccs_int     index_p      - the cell index.
    !> @param[out] cell_locator loc_p - the cell locator object linking a cell index with
    !!                                          the mesh.
    module subroutine set_cell_location(mesh, index_p, loc_p)
      type(ccs_mesh), target, intent(in) :: mesh
      integer(ccs_int), intent(in) :: index_p
      type(cell_locator), intent(out) :: loc_p
    end subroutine set_cell_location

    !>  Constructs a face locator object.
    !
    !> @description Creates the association between a face relative to a cell, i.e. to access the
    !!              nth face of cell i.
    !
    !> @param[in]  mesh         mesh      - the mesh object being referred to.
    !> @param[in]  ccs_int     index_p      - the index of the cell whose face is being accessed.
    !> @param[in]  ccs_int     cell_face_ctr - the cell-local index of the face.
    !> @param[out] face_locator loc_f - the face locator object linking a cell-relative
    !!                                          index with the mesh.
    module subroutine set_face_location(mesh, index_p, cell_face_ctr, loc_f)
      type(ccs_mesh), target, intent(in) :: mesh
      integer(ccs_int), intent(in) :: index_p
      integer(ccs_int), intent(in) :: cell_face_ctr
      type(face_locator), intent(out) :: loc_f
    end subroutine set_face_location

    !>  Constructs a neighbour locator object.
    !
    !> @description Creates the association between a neighbour cell F relative to cell P, i.e. to
    !!              access the nth neighbour of cell i.
    !
    !> @param[in]  cell_locator      loc_p      - the cell locator object of the cell whose
    !!                                                    neighbour is being accessed.
    !> @param[in]  ccs_int           nb_counter - the cell-local index of the neighbour.
    !> @param[out] neighbour_locator loc_nb - the neighbour locator object linking a
    !!                                                    cell-relative index with the mesh.
    module subroutine set_neighbour_location(loc_p, nb_counter, loc_nb)
      type(cell_locator), intent(in) :: loc_p
      integer(ccs_int), intent(in) :: nb_counter
      type(neighbour_locator), intent(out) :: loc_nb
    end subroutine set_neighbour_location

    !>  Set face index
    module subroutine set_face_index(index_p, cell_face_ctr, index_f, mesh)
      integer(ccs_int), intent(in) :: index_p
      integer(ccs_int), intent(in) :: cell_face_ctr
      integer(ccs_int), intent(in) :: index_f
      type(ccs_mesh), target, intent(inout) :: mesh
    end subroutine set_face_index

    !>  Returns the normal vector of a face
    !
    !> @param[in]  face_locator    loc_f - the face locator object.
    !> @param[out] real(ccs_real) normal(ndim)  - an ndimensional array representing the face normal
    !!                                             vector.
    module subroutine get_face_normal(loc_f, normal)
      type(face_locator), intent(in) :: loc_f
      real(ccs_real), dimension(ndim), intent(out) :: normal
    end subroutine get_face_normal

    !>  Returns the area of a face
    !
    !> @param[in]  face_locator    loc_f - the face locator object.
    !> @param[out] real(ccs_real) area          - the face area.
    module subroutine get_face_area(loc_f, area)
      type(face_locator), intent(in) :: loc_f
      real(ccs_real), intent(out) :: area
    end subroutine get_face_area

    !>  Returns the centre of a cell
    !
    !> @param[in]  cell_locator     loc_p - the cell locator object.
    !> @param[out] real(ccs_real)  x(ndim)       - an ndimensional array representing the cell centre.
    module subroutine get_cell_centre(loc_p, x)
      type(cell_locator), intent(in) :: loc_p
      real(ccs_real), dimension(ndim), intent(out) :: x
    end subroutine get_cell_centre

    !>  Returns the centre of a neighbour cell
    !
    !> @param[in]  neighbour_locator loc_nb - the neighbour locator object.
    !> @param[out] real(ccs_real)   x(ndim)            - an ndimensional array representing the
    !!                                                    neighbour cell centre.
    module subroutine get_neighbour_centre(loc_nb, x)
      type(neighbour_locator), intent(in) :: loc_nb
      real(ccs_real), dimension(ndim), intent(out) :: x
    end subroutine get_neighbour_centre

    !>  Returns the centre of a face
    !
    !> @param[in]  face_locator     loc_f - the face locator object.
    !> @param[out] real(ccs_real)  x(ndim)       - an ndimensional array representing the face centre.
    module subroutine get_face_centre(loc_f, x)
      type(face_locator), intent(in) :: loc_f
      real(ccs_real), dimension(ndim), intent(out) :: x
    end subroutine get_face_centre

    !>  Returns the volume of a cell
    !
    !> @param[in] cell_locator     loc_p - the cell locator object.
    !> @param[out] real(ccs_real) V             - the cell volume.
    module subroutine get_cell_volume(loc_p, V)
      type(cell_locator), intent(in) :: loc_p
      real(ccs_real), intent(out) :: V
    end subroutine get_cell_volume

    !>  Returns the volume of a neighbour cell
    !
    !> @param[in] neighbour_locator loc_nb - the neighbour locator object.
    !> @param[out] real(ccs_real)  V                  - the neighbour cell volume.
    module subroutine get_neighbour_volume(loc_nb, V)
      type(neighbour_locator), intent(in) :: loc_nb
      real(ccs_real), intent(out) :: V
    end subroutine get_neighbour_volume

    !>  Returns the global index of a cell
    !
    !> @param[in]  cell_locator      loc_p - the cell locator object.
    !> @param[out] integer(ccs_int) global_index_p          - the global index of the cell.
    module subroutine get_cell_global_index(loc_p, global_index_p)
      type(cell_locator), intent(in) :: loc_p
      integer(ccs_int), intent(out) :: global_index_p
    end subroutine get_cell_global_index

    !>  Returns the global index of a neighbouring cell
    !
    !> @param[in]  neighbour_locator loc_nb - the neighbour locator object.
    !> @param[out] integer(ccs_int) global_index_nb             - the global index of the neighbour cell.
    module subroutine get_neighbour_global_index(loc_nb, global_index_nb)
      type(neighbour_locator), intent(in) :: loc_nb
      integer(ccs_int), intent(out) :: global_index_nb
    end subroutine get_neighbour_global_index

    !>  Returns the neighbour count of a cell (including boundary neighbours)
    !
    !> @param[in]  cell_locator      loc_p - the cell locator object.
    !> @param[out] integer(ccs_int) nnb           - the neighbour count of the cell.
    module subroutine cell_count_neighbours(loc_p, nnb)
      type(cell_locator), intent(in) :: loc_p
      integer(ccs_int), intent(out) :: nnb
    end subroutine cell_count_neighbours

    !>  Returns the boundary status of a neighbouring cell
    !
    !> @param[in]  neighbour_locator loc_nb - the neighbour locator object.
    !> @param[out] logical           is_boundary        - the boundary status of the neighbour.
    module subroutine get_neighbour_boundary_status(loc_nb, is_boundary)
      type(neighbour_locator), intent(in) :: loc_nb
      logical, intent(out) :: is_boundary
    end subroutine get_neighbour_boundary_status

    !>  Returns the boundary status of a face
    !
    !> @param[in]  face_locator loc_f - the face locator object.
    !> @param[out] logical      is_boundary   - the boundary status of the neighbour.
    module subroutine get_face_boundary_status(loc_f, is_boundary)
      type(face_locator), intent(in) :: loc_f
      logical, intent(out) :: is_boundary
    end subroutine get_face_boundary_status

    !>  Returns the local distribution status of a neighbouring cell
    !
    !> @description Given a distributed mesh, a processor needs both the cells within its partition
    !!              and cells from the surrounding halo - this subroutine get_indicates whether a
    !!              cell's neighbour is within the local partition or the halo.
    !
    !> @param[in]  neighbour_locator loc_nb - the neighbour locator object.
    !> @param[out] logical           is_local           - the local status of the neighbour.
    module subroutine get_local_status(loc_nb, is_local)
      type(neighbour_locator), intent(in) :: loc_nb
      logical, intent(out) :: is_local
    end subroutine get_local_status

    !>  Returns the local index of a cell
    !
    !> @description Generally the local index of a cell is should be the same as its location within
    !!              the local cell vector - this particular subroutine get_is therefore expected of
    !!              limited use and is mostly present for uniformity.
    !
    !> @param[in]  cell_locator      loc_p - the cell locator object.
    !> @param[out] integer(ccs_int) index_p           - the local index of the cell.
    module subroutine get_cell_local_index(loc_p, index_p)
      type(cell_locator), intent(in) :: loc_p
      integer(ccs_int), intent(out) :: index_p
    end subroutine get_cell_local_index

    !>  Returns the local index of a neighbouring cell
    !
    !> @param[in]  neighbour_locator loc_nb - the neighbour locator object.
    !> @param[out] integer(ccs_int) index_nb              - the local index of the neighbour cell.
    module subroutine get_neighbour_local_index(loc_nb, index_nb)
      type(neighbour_locator), intent(in) :: loc_nb
      integer(ccs_int), intent(out) :: index_nb
    end subroutine get_neighbour_local_index

    !>  Returns the local index of a face
    !
    !> @param[in]  face_locator     loc_f     - the face locator object.
    !> @param[out] integer(ccs_int) index_f   - the local index of the face.
    module subroutine get_face_local_index(loc_f, index_f)
      type(face_locator), intent(in) :: loc_f
      integer(ccs_int), intent(out) :: index_f
    end subroutine get_face_local_index

    !>  Returns the distance between two cell centres
    !
    !> @param[in]  cell_locator      loc_p  - The cell distance is measured from.
    !> @param[in]  neighbour_locator loc_nb - The cell distance is measured to.
    !> @param[out] ccs_real         dx     - ndim-array of the distance
    module subroutine get_neighbour_distance(loc_p, loc_nb, dx)
      type(cell_locator), intent(in) :: loc_p
      type(neighbour_locator), intent(in) :: loc_nb
      real(ccs_real), dimension(ndim), intent(out) :: dx
    end subroutine get_neighbour_distance

    !>  Returns the distance from cell to face centres
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

