!v Module file meshing.mod
!
!  Module defining meshing interface for ASiMoV-CCS

module meshing

  use constants, only: ndim
  use kinds, only: ccs_int, ccs_real
  use types, only: ccs_mesh, face_locator, cell_locator, neighbour_locator, vert_locator

  implicit none

  private
  public :: set_face_location
  public :: set_cell_location
  public :: set_neighbour_location
  public :: set_vert_location
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
  public :: set_centre
  public :: set_area
  
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

  interface set_centre
     module procedure set_cell_centre
     module procedure set_face_centre
     module procedure set_vert_centre
  end interface set_centre
  
  interface

    !v Constructs a cell locator object.
    !
    !  Creates the association between a mesh and cell index, storing it in the
    !  returned cell locator object.
    module subroutine set_cell_location(mesh, index_p, loc_p)
      type(ccs_mesh), target, intent(in) :: mesh !< the mesh object being referred to.
      integer(ccs_int), intent(in) :: index_p    !< the cell index.
      type(cell_locator), intent(out) :: loc_p   !< the cell locator object linking a cell index with teh mesh.
    end subroutine set_cell_location

    !v Constructs a face locator object.
    !
    !  Creates the association between a face relative to a cell, i.e. to access the
    !  nth face of cell i.
    module subroutine set_face_location(mesh, index_p, cell_face_ctr, loc_f)
      type(ccs_mesh), target, intent(in) :: mesh    !< the mesh object being referred to.
      integer(ccs_int), intent(in) :: index_p       !< the index of the cell whose face is being accessed.
      integer(ccs_int), intent(in) :: cell_face_ctr !< the cell-local index of the face.
      type(face_locator), intent(out) :: loc_f      !< the face locator object linking a cell-relative index with the mesh.
    end subroutine set_face_location

    !v Constructs a neighbour locator object.
    !
    !  Creates the association between a neighbour cell F relative to cell P, i.e. to
    !  access the nth neighbour of cell i.
    module subroutine set_neighbour_location(loc_p, nb_counter, loc_nb)
      type(cell_locator), intent(in) :: loc_p        !< the cell locator object of the cell whose neighbour is being accessed.
      integer(ccs_int), intent(in) :: nb_counter     !< the cell-local index of the neighbour.
      type(neighbour_locator), intent(out) :: loc_nb !< the neighbour locator object linking a cell-relative index with the mesh.
    end subroutine set_neighbour_location

    !v Constructs a vertex locator object.
    !
    !  Creates the association between a vertex relative to a cell, i.e. to access the
    !  nth vertex of cell i.
    module subroutine set_vert_location(mesh, index_p, cell_vert_ctr, loc_v)
      type(ccs_mesh), target, intent(in) :: mesh    !< the mesh object being referred to.
      integer(ccs_int), intent(in) :: index_p       !< the index of the cell whose vertex is being accessed.
      integer(ccs_int), intent(in) :: cell_vert_ctr !< the cell-local index of the vertex.
      type(vert_locator), intent(out) :: loc_v      !< the vertex locator object linking a cell-relative index with the mesh.
    end subroutine set_vert_location

    !> Set face index
    module subroutine set_face_index(index_p, cell_face_ctr, index_f, mesh)
      integer(ccs_int), intent(in) :: index_p
      integer(ccs_int), intent(in) :: cell_face_ctr
      integer(ccs_int), intent(in) :: index_f
      type(ccs_mesh), target, intent(inout) :: mesh
    end subroutine set_face_index

    !> Returns the normal vector of a face
    module subroutine get_face_normal(loc_f, normal)
      type(face_locator), intent(in) :: loc_f                !< the face locator object.
      real(ccs_real), dimension(ndim), intent(out) :: normal !< an ndimensional array representing the face normal vector.
    end subroutine get_face_normal

    !> Returns the area of a face
    module subroutine get_face_area(loc_f, area)
      type(face_locator), intent(in) :: loc_f !< the face locator object.
      real(ccs_real), intent(out) :: area     !< the face area.
    end subroutine get_face_area

    !> Returns the centre of a cell
    module subroutine get_cell_centre(loc_p, x)
      type(cell_locator), intent(in) :: loc_p           !< the cell locator object.
      real(ccs_real), dimension(:), intent(out) :: x !< an ndimensional array representing the cell centre.
    end subroutine get_cell_centre
    
    !> Returns the centre of a neighbour cell
    module subroutine get_neighbour_centre(loc_nb, x)
      type(neighbour_locator), intent(in) :: loc_nb     !< the neighbour locator object.
      real(ccs_real), dimension(ndim), intent(out) :: x !< an ndimensional array representing the neighbour cell centre.
    end subroutine get_neighbour_centre

    !> Returns the centre of a face
    module subroutine get_face_centre(loc_f, x)
      type(face_locator), intent(in) :: loc_f           !< the face locator object.
      real(ccs_real), dimension(ndim), intent(out) :: x !< an ndimensional array representing the face centre.
    end subroutine get_face_centre

    !> Returns the volume of a cell
    module subroutine get_cell_volume(loc_p, V)
      type(cell_locator), intent(in) :: loc_p !< the cell locator object.
      real(ccs_real), intent(out) :: V        !< the cell volume.
    end subroutine get_cell_volume

    !> Returns the volume of a neighbour cell
    module subroutine get_neighbour_volume(loc_nb, V)
      type(neighbour_locator), intent(in) :: loc_nb !< the neighbour locator object.
      real(ccs_real), intent(out) :: V              !< the neighbour cell volume.
    end subroutine get_neighbour_volume

    !> Returns the global index of a cell
    module subroutine get_cell_global_index(loc_p, global_index_p)
      type(cell_locator), intent(in) :: loc_p         !< the cell locator object.
      integer(ccs_int), intent(out) :: global_index_p !< the global index of the cell.
    end subroutine get_cell_global_index

    !> Returns the global index of a neighbouring cell
    module subroutine get_neighbour_global_index(loc_nb, global_index_nb)
      type(neighbour_locator), intent(in) :: loc_nb    !< the neighbour locator object.
      integer(ccs_int), intent(out) :: global_index_nb !< the global index of the neighbour cell.
    end subroutine get_neighbour_global_index

    !> Returns the neighbour count of a cell (including boundary neighbours)
    module subroutine cell_count_neighbours(loc_p, nnb)
      type(cell_locator), intent(in) :: loc_p !< the cell locator object.
      integer(ccs_int), intent(out) :: nnb    !< the neighbour count of the cell.
    end subroutine cell_count_neighbours

    !> Returns the boundary status of a neighbouring cell
    module subroutine get_neighbour_boundary_status(loc_nb, is_boundary)
      type(neighbour_locator), intent(in) :: loc_nb !< the neighbour locator object.
      logical, intent(out) :: is_boundary           !< the boundary status of the neighbour.
    end subroutine get_neighbour_boundary_status

    !> Returns the boundary status of a face
    module subroutine get_face_boundary_status(loc_f, is_boundary)
      type(face_locator), intent(in) :: loc_f !< the face locator object.
      logical, intent(out) :: is_boundary     !< the boundary status of the neighbour.
    end subroutine get_face_boundary_status

    !v Returns the local distribution status of a neighbouring cell
    !
    !  Given a distributed mesh, a processor needs both the cells within its partition
    !  and cells from the surrounding halo - this subroutine get_indicates whether a
    !  cell's neighbour is within the local partition or the halo.
    module subroutine get_local_status(loc_nb, is_local)
      type(neighbour_locator), intent(in) :: loc_nb !< the neighbour locator object.
      logical, intent(out) :: is_local !< the local status of the neighbour.
    end subroutine get_local_status

    !v Returns the local index of a cell
    !
    !  Generally the local index of a cell is should be the same as its location within
    !  the local cell vector - this particular subroutine get_is therefore expected of
    !  limited use and is mostly present for uniformity.
    module subroutine get_cell_local_index(loc_p, index_p)
      type(cell_locator), intent(in) :: loc_p  !< the cell locator object.
      integer(ccs_int), intent(out) :: index_p !< the local index of the cell.
    end subroutine get_cell_local_index

    !> Returns the local index of a neighbouring cell
    module subroutine get_neighbour_local_index(loc_nb, index_nb)
      type(neighbour_locator), intent(in) :: loc_nb !< the neighbour locator object.
      integer(ccs_int), intent(out) :: index_nb     !< the local index of the neighbour cell.
    end subroutine get_neighbour_local_index

    !> Returns the local index of a face
    module subroutine get_face_local_index(loc_f, index_f)
      type(face_locator), intent(in) :: loc_f !< the face locator object.
      integer(ccs_int), intent(out) :: index_f !< the local index of the face.
    end subroutine get_face_local_index

    !> Returns the distance between two cell centres
    module subroutine get_neighbour_distance(loc_p, loc_nb, dx)
      type(cell_locator), intent(in) :: loc_p            !< The cell distance is measured from.
      type(neighbour_locator), intent(in) :: loc_nb      !< The cell distance is measured to.
      real(ccs_real), dimension(ndim), intent(out) :: dx !< ndim-array of the distance
    end subroutine get_neighbour_distance

    !> Returns the distance from cell to face centres
    module subroutine get_face_distance(loc_p, loc_f, dx)
      type(cell_locator), intent(in) :: loc_p            !< The cell distance is measured from.
      type(face_locator), intent(in) :: loc_f            !< The face distance is measured to.
      real(ccs_real), dimension(ndim), intent(out) :: dx !< ndim-array of the distance
    end subroutine get_face_distance

    !> Set the cell centre of specified cell
    module subroutine set_cell_centre(loc_p, x_p)
      type(cell_locator), intent(in) :: loc_p         !< The cell locator object.
      real(ccs_real), dimension(:), intent(in) :: x_p !< The cell centre array.
    end subroutine set_cell_centre

    !> Set the face centre of specified face
    module subroutine set_face_centre(loc_f, x_f)
      type(face_locator), intent(in) :: loc_f         !< The face locator object.
      real(ccs_real), dimension(:), intent(in) :: x_f !< The face centre array.
    end subroutine set_face_centre

    !> Set the centre of specified vertex
    module subroutine set_vert_centre(loc_v, x_v)
      type(vert_locator), intent(in) :: loc_v         !< The vertex locator object.
      real(ccs_real), dimension(:), intent(in) :: x_v !< The vertex centre array.
    end subroutine set_vert_centre

    !> Set the area of specified face
    module subroutine set_area(area, loc_f)
      real(ccs_real), intent(in) :: area      !< The face area
      type(face_locator), intent(in) :: loc_f !< The face locator object
    end subroutine set_area
    
  end interface

end module meshing

