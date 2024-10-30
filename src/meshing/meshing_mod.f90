!v Module file meshing.mod
!
!  Module defining meshing interface for ASiMoV-CCS

module meshing

  use constants, only: ndim
  use kinds, only: ccs_int, ccs_real, ccs_long
  use types, only: ccs_mesh, topology, face_locator, &
                   cell_locator, neighbour_locator, vert_locator

  implicit none

  private
  public :: set_mesh_object
  public :: nullify_mesh_object
  public :: is_mesh_set
  public :: set_topo_object
  public :: nullify_topo_object
  public :: is_topo_set
  public :: create_face_locator
  public :: create_cell_locator
  public :: create_neighbour_locator
  public :: create_vert_locator
  public :: set_face_index
  public :: get_face_normal
  public :: get_face_area
  public :: get_centre
  public :: get_volume
  public :: get_global_index, set_global_index
  public :: get_natural_index, set_natural_index
  public :: get_local_index, set_local_index
  public :: get_boundary_status
  public :: get_local_status
  public :: count_neighbours
  public :: get_distance
  public :: get_local_num_cells, set_local_num_cells
  public :: get_total_num_cells, set_total_num_cells
  public :: get_global_num_cells, set_global_num_cells
  public :: get_halo_num_cells, set_halo_num_cells
  public :: get_global_num_faces, set_global_num_faces
  public :: get_num_faces, set_num_faces
  public :: get_max_faces, set_max_faces
  public :: get_global_num_vertices, set_global_num_vertices
  public :: get_vert_per_cell, set_vert_per_cell
  public :: set_centre
  public :: set_area
  public :: set_volume
  public :: set_normal
  public :: get_face_interpolation
  public :: set_face_interpolation
  public :: get_mesh_generated, set_mesh_generated
  public :: get_bc_id
  public :: find_entities
  
  interface get_centre
    module procedure get_cell_centre
    module procedure get_neighbour_centre
    module procedure get_face_centre
    module procedure get_vert_centre
  end interface get_centre

  interface get_global_index
    module procedure get_cell_global_index
    module procedure get_neighbour_global_index
  end interface get_global_index

  interface get_natural_index
    module procedure get_cell_natural_index
    module procedure get_neighbour_natural_index
  end interface get_natural_index

  interface set_global_index
    module procedure set_face_global_index
    module procedure set_cell_global_index
  end interface set_global_index

  interface set_natural_index
    module procedure set_cell_natural_index
  end interface set_natural_index

  interface get_local_index
    module procedure get_cell_local_index
    module procedure get_neighbour_local_index
    module procedure get_face_local_index
  end interface get_local_index

  interface set_local_index
    module procedure set_neighbour_local_index
  end interface set_local_index

  interface count_neighbours
    module procedure get_cell_count_neighbours
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
    module procedure get_face_neighbour_distance
  end interface get_distance

  interface get_local_num_cells
    module procedure get_local_num_cells_int
    module procedure get_local_num_cells_long
  end interface get_local_num_cells

  interface set_centre
    module procedure set_cell_centre
    module procedure set_face_centre
    module procedure set_vert_centre
  end interface set_centre

  interface create_neighbour_locator
    module procedure create_face_neighbour_locator
  end interface create_neighbour_locator

  interface get_local_status
    module procedure get_neighbour_local_status
  end interface get_local_status

  interface find_entities
    module procedure find_face_entities
  end interface find_entities
  
  interface

    module subroutine set_mesh_object(input_mesh)
      type(ccs_mesh), target, intent(inout) :: input_mesh           !< The mesh
    end subroutine

    module subroutine nullify_mesh_object()
    end subroutine

    module subroutine set_topo_object(input_topo)
      type(topology), target, intent(inout) :: input_topo           !< The mesh topology object
    end subroutine

    module subroutine nullify_topo_object()
    end subroutine

    pure module function is_mesh_set()
      logical :: is_mesh_set
    end function

    pure module function is_topo_set()
      logical :: is_topo_set
    end function

    !v Constructs a cell locator object.
    !
    !  Creates the association between a mesh and cell index, storing it in the
    !  returned cell locator object.
    pure module subroutine create_cell_locator(index_p, loc_p)
      integer(ccs_int), intent(in) :: index_p    !< the cell index.
      type(cell_locator), intent(out) :: loc_p   !< the cell locator object linking a cell index with teh mesh.
    end subroutine create_cell_locator

    !v Constructs a face locator object.
    !
    !  Creates the association between a face relative to a cell, i.e. to access the
    !  nth face of cell i.
    pure module subroutine create_face_locator(index_p, cell_face_ctr, loc_f)
      integer(ccs_int), intent(in) :: index_p       !< the index of the cell whose face is being accessed.
      integer(ccs_int), intent(in) :: cell_face_ctr !< the cell-local index of the face.
      type(face_locator), intent(out) :: loc_f      !< the face locator object linking a cell-relative index with the mesh.
    end subroutine create_face_locator

    !v Constructs a neighbour locator object.
    !
    !  Creates the association between a neighbour cell F relative to cell P, i.e. to
    !  access the nth neighbour of cell i.
    pure module subroutine create_face_neighbour_locator(loc_p, nb_counter, loc_nb)
      type(cell_locator), intent(in) :: loc_p        !< the cell locator object of the cell whose neighbour is being accessed.
      integer(ccs_int), intent(in) :: nb_counter     !< the cell-local index of the neighbour.
      type(neighbour_locator), intent(out) :: loc_nb !< the neighbour locator object linking a cell-relative index with the mesh.
    end subroutine create_face_neighbour_locator

    !v Constructs a vertex locator object.
    !
    !  Creates the association between a vertex relative to a cell, i.e. to access the
    !  nth vertex of cell i.
    pure module subroutine create_vert_locator(index_p, cell_vert_ctr, loc_v)
      integer(ccs_int), intent(in) :: index_p       !< the index of the cell whose vertex is being accessed.
      integer(ccs_int), intent(in) :: cell_vert_ctr !< the cell-local index of the vertex.
      type(vert_locator), intent(out) :: loc_v      !< the vertex locator object linking a cell-relative index with the mesh.
    end subroutine create_vert_locator

    !> Set face index
    module subroutine set_face_index(index_p, cell_face_ctr, index_f)
      integer(ccs_int), intent(in) :: index_p
      integer(ccs_int), intent(in) :: cell_face_ctr
      integer(ccs_int), intent(in) :: index_f
    end subroutine set_face_index

    !> Returns the normal vector of a face
    pure module subroutine get_face_normal(loc_f, normal)
      type(face_locator), intent(in) :: loc_f                !< the face locator object.
      real(ccs_real), dimension(ndim), intent(out) :: normal !< an ndimensional array representing the face normal vector.
    end subroutine get_face_normal

    !> Returns the area of a face
    pure module subroutine get_face_area(loc_f, area)
      type(face_locator), intent(in) :: loc_f !< the face locator object.
      real(ccs_real), intent(out) :: area     !< the face area.
    end subroutine get_face_area

    !> Returns the centre of a cell
    pure module subroutine get_cell_centre(loc_p, x)
      type(cell_locator), intent(in) :: loc_p           !< the cell locator object.
      real(ccs_real), dimension(:), intent(out) :: x !< an ndimensional array representing the cell centre.
    end subroutine get_cell_centre

    !> Returns the centre of a neighbour cell
    pure module subroutine get_neighbour_centre(loc_nb, x)
      type(neighbour_locator), intent(in) :: loc_nb     !< the neighbour locator object.
      real(ccs_real), dimension(ndim), intent(out) :: x !< an ndimensional array representing the neighbour cell centre.
    end subroutine get_neighbour_centre

    !> Returns the centre of a face
    pure module subroutine get_face_centre(loc_f, x)
      type(face_locator), intent(in) :: loc_f           !< the face locator object.
      real(ccs_real), dimension(ndim), intent(out) :: x !< an ndimensional array representing the face centre.
    end subroutine get_face_centre

    !> Returns the centre of a vertex
    pure module subroutine get_vert_centre(loc_v, x)
      type(vert_locator), intent(in) :: loc_v           !< the vertex locator object.
      real(ccs_real), dimension(:), intent(out) :: x !< an ndimensional array representing the vertex centre.
    end subroutine get_vert_centre

    !> Returns the volume of a cell
    pure module subroutine get_cell_volume(loc_p, V)
      type(cell_locator), intent(in) :: loc_p !< the cell locator object.
      real(ccs_real), intent(out) :: V        !< the cell volume.
    end subroutine get_cell_volume

    !> Returns the volume of a neighbour cell
    pure module subroutine get_neighbour_volume(loc_nb, V)
      type(neighbour_locator), intent(in) :: loc_nb !< the neighbour locator object.
      real(ccs_real), intent(out) :: V              !< the neighbour cell volume.
    end subroutine get_neighbour_volume

    !> Returns the global index of a cell
    pure module subroutine get_cell_global_index(loc_p, global_index_p)
      type(cell_locator), intent(in) :: loc_p         !< the cell locator object.
      integer(ccs_int), intent(out) :: global_index_p !< the global index of the cell.
    end subroutine get_cell_global_index

    !> Sets the global index of a cell
    module subroutine set_cell_global_index(global_index_p, loc_p)
      integer(ccs_int), intent(in) :: global_index_p !< the global index of the cell.
      type(cell_locator), intent(inout) :: loc_p     !< the cell locator object.
    end subroutine set_cell_global_index

    !v Returns the natural index of a cell
    !
    ! @note@ The natural index is the original global index, whereas the global index indicates the
    !        indexing in the current ordering.
    pure module subroutine get_cell_natural_index(loc_p, natural_index_p)
      type(cell_locator), intent(in) :: loc_p          !< the cell locator object.
      integer(ccs_int), intent(out) :: natural_index_p !< the natural index of the cell.
    end subroutine get_cell_natural_index

    !v Sets the natural index of a cell
    !
    ! @note@ The natural index is the original global index, whereas the global index indicates the
    !        indexing in the current ordering.
    module subroutine set_cell_natural_index(natural_index_p, loc_p)
      integer(ccs_int), intent(in) :: natural_index_p !< the natural index of the cell.
      type(cell_locator), intent(inout) :: loc_p      !< the cell locator object.
    end subroutine set_cell_natural_index

    !> Returns the global index of a neighbouring cell
    pure module subroutine get_neighbour_global_index(loc_nb, global_index_nb)
      type(neighbour_locator), intent(in) :: loc_nb    !< the neighbour locator object.
      integer(ccs_int), intent(out) :: global_index_nb !< the global index of the neighbour cell.
    end subroutine get_neighbour_global_index

    !> Returns the natural index of a neighbouring cell
    pure module subroutine get_neighbour_natural_index(loc_nb, natural_index_nb)
      type(neighbour_locator), intent(in) :: loc_nb     !< the neighbour locator object.
      integer(ccs_int), intent(out) :: natural_index_nb !< the natural index of the neighbour cell.
    end subroutine get_neighbour_natural_index

    !> Sets the global index of a face
    module subroutine set_face_global_index(global_index_f, loc_f)
      integer(ccs_int), intent(in) :: global_index_f !< The global index of the face.
      type(face_locator), intent(inout) :: loc_f     !< The face locator object.
    end subroutine set_face_global_index

    !> Returns the neighbour count of a cell (including boundary neighbours)
    pure module subroutine get_cell_count_neighbours(loc_p, nnb)
      type(cell_locator), intent(in) :: loc_p !< the cell locator object.
      integer(ccs_int), intent(out) :: nnb    !< the neighbour count of the cell.
    end subroutine get_cell_count_neighbours

    !> Returns the boundary status of a neighbouring cell
    pure module subroutine get_neighbour_boundary_status(loc_nb, is_boundary)
      type(neighbour_locator), intent(in) :: loc_nb !< the neighbour locator object.
      logical, intent(out) :: is_boundary           !< the boundary status of the neighbour.
    end subroutine get_neighbour_boundary_status

    !> Returns the boundary status of a face
    pure module subroutine get_face_boundary_status(loc_f, is_boundary)
      type(face_locator), intent(in) :: loc_f !< the face locator object.
      logical, intent(out) :: is_boundary     !< the boundary status of the neighbour.
    end subroutine get_face_boundary_status

    !v Returns the local distribution status of a neighbouring cell
    !
    !  Given a distributed mesh, a processor needs both the cells within its partition
    !  and cells from the surrounding halo - this subroutine get_indicates whether a
    !  cell's neighbour is within the local partition or the halo.
    pure module subroutine get_neighbour_local_status(loc_nb, is_local)
      type(neighbour_locator), intent(in) :: loc_nb !< the neighbour locator object.
      logical, intent(out) :: is_local !< the local status of the neighbour.
    end subroutine get_neighbour_local_status

    !v Returns the local index of a cell
    !
    !  Generally the local index of a cell is should be the same as its location within
    !  the local cell vector - this particular subroutine get_is therefore expected of
    !  limited use and is mostly present for uniformity.
    pure module subroutine get_cell_local_index(loc_p, index_p)
      type(cell_locator), intent(in) :: loc_p  !< the cell locator object.
      integer(ccs_int), intent(out) :: index_p !< the local index of the cell.
    end subroutine get_cell_local_index

    !> Returns the local index of a neighbouring cell
    pure module subroutine get_neighbour_local_index(loc_nb, index_nb)
      type(neighbour_locator), intent(in) :: loc_nb !< the neighbour locator object.
      integer(ccs_int), intent(out) :: index_nb     !< the local index of the neighbour cell.
    end subroutine get_neighbour_local_index

    !> Sets the local index of a neighbouring cell
    module subroutine set_neighbour_local_index(index_nb, loc_nb)
      integer(ccs_int), intent(in) :: index_nb     !< the local index of the neighbour cell.
      type(neighbour_locator), intent(inout) :: loc_nb !< the neighbour locator object.
    end subroutine set_neighbour_local_index

    !> Returns the local index of a face
    pure module subroutine get_face_local_index(loc_f, index_f)
      type(face_locator), intent(in) :: loc_f !< the face locator object.
      integer(ccs_int), intent(out) :: index_f !< the local index of the face.
    end subroutine get_face_local_index

    !> Returns the distance between two cell centres
    pure module subroutine get_neighbour_distance(loc_p, loc_nb, dx)
      type(cell_locator), intent(in) :: loc_p            !< The cell distance is measured from.
      type(neighbour_locator), intent(in) :: loc_nb      !< The cell distance is measured to.
      real(ccs_real), dimension(ndim), intent(out) :: dx !< ndim-array of the distance
    end subroutine get_neighbour_distance

    !> Returns the distance from cell to face centres
    pure module subroutine get_face_distance(loc_p, loc_f, dx)
      type(cell_locator), intent(in) :: loc_p            !< The cell distance is measured from.
      type(face_locator), intent(in) :: loc_f            !< The face distance is measured to.
      real(ccs_real), dimension(ndim), intent(out) :: dx !< ndim-array of the distance
    end subroutine get_face_distance

    !> Returns the distance from neighbour cell to face centres
    pure module subroutine get_face_neighbour_distance(loc_nb, loc_f, dx)
      type(neighbour_locator), intent(in) :: loc_nb      !< The cell distance is measured from.
      type(face_locator), intent(in) :: loc_f            !< The face distance is measured to.
      real(ccs_real), dimension(ndim), intent(out) :: dx !< ndim-array of the distance
    end subroutine get_face_neighbour_distance

    !> Sets the mesh topology local cell count
    module subroutine set_local_num_cells(local_num_cells)
      integer(ccs_int), intent(in) :: local_num_cells !< The local cell count
    end subroutine set_local_num_cells

    !> Gets the mesh topology local cell count.
    pure module subroutine get_local_num_cells_int(local_num_cells)
      integer(ccs_int), intent(out) :: local_num_cells !< The local cell count
    end subroutine get_local_num_cells_int
    !v Gets the mesh topology local cell count.
    !
    !  Handles case when using a long integer to access the internal topology data.
    pure module subroutine get_local_num_cells_long(local_num_cells)
      integer(ccs_long), intent(out) :: local_num_cells !< The local cell count
    end subroutine get_local_num_cells_long

    !> Sets the mesh topology total cell count.
    module subroutine set_total_num_cells(total_num_cells)
      integer(ccs_int), intent(in) :: total_num_cells !< The total cell count
    end subroutine set_total_num_cells

    !> Gets the mesh total cell count.
    pure module subroutine get_total_num_cells(total_num_cells)
      integer(ccs_int), intent(out) :: total_num_cells !< The total cell count
    end subroutine get_total_num_cells

    !> Sets the mesh global cell count.
    module subroutine set_global_num_cells(global_num_cells)
      integer(ccs_int), intent(in) :: global_num_cells !< The global cell count
    end subroutine set_global_num_cells

    !> Gets the mesh topology global cell count.
    pure module subroutine get_global_num_cells(global_num_cells)
      integer(ccs_int), intent(out) :: global_num_cells !< The global cell count
    end subroutine get_global_num_cells

    !> Sets the mesh halo cell count.
    module subroutine set_halo_num_cells(halo_num_cells)
      integer(ccs_int), intent(in) :: halo_num_cells !< The halo cell count
    end subroutine set_halo_num_cells

    !> Gets the mesh halo cell count.
    pure module subroutine get_halo_num_cells(halo_num_cells)
      integer(ccs_int), intent(out) :: halo_num_cells !< The halo cell count
    end subroutine get_halo_num_cells

    !> Sets the mesh global face count.
    module subroutine set_global_num_faces(global_num_faces)
      integer(ccs_int), intent(in) :: global_num_faces !< The global face count
    end subroutine set_global_num_faces

    !> Gets the mesh topology global face count.
    pure module subroutine get_global_num_faces(global_num_faces)
      integer(ccs_int), intent(out) :: global_num_faces !< The global face count
    end subroutine get_global_num_faces

    !> Sets the mesh face count.
    module subroutine set_num_faces(num_faces)
      integer(ccs_int), intent(in) :: num_faces !< The face count
    end subroutine set_num_faces

    !> Gets the mesh face count.
    pure module subroutine get_num_faces(num_faces)
      integer(ccs_int), intent(out) :: num_faces !< The face count
    end subroutine get_num_faces

    !> Sets the mesh face count.
    module subroutine set_max_faces(max_faces)
      integer(ccs_int), intent(in) :: max_faces !< The face count
    end subroutine set_max_faces

    !> Gets the mesh topology face count.
    pure module subroutine get_max_faces(max_faces)
      integer(ccs_int), intent(out) :: max_faces !< The face count
    end subroutine get_max_faces

    !> Sets the global number of vertices.
    module subroutine set_global_num_vertices(global_num_vertices)
      integer(ccs_int), intent(in) :: global_num_vertices !< The global number of vertices
    end subroutine set_global_num_vertices

    !> Gets the global number of vertices.
    pure module subroutine get_global_num_vertices(global_num_vertices)
      integer(ccs_int), intent(out) :: global_num_vertices !< The global number of vertices
    end subroutine get_global_num_vertices

    !> Sets the number of vertices per cell.
    module subroutine set_vert_per_cell(vert_per_cell)
      integer(ccs_int), intent(in) :: vert_per_cell !< The number of vertices per cell
    end subroutine set_vert_per_cell

    !> Gets the number of vertices per cell.
    pure module subroutine get_vert_per_cell(vert_per_cell)
      integer(ccs_int), intent(out) :: vert_per_cell !< The number of vertices per cell
    end subroutine get_vert_per_cell
    
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

    !> Set the volume of specified cell
    module subroutine set_volume(volume, loc_p)
      real(ccs_real), intent(in) :: volume      !< The cell volume
      type(cell_locator), intent(in) :: loc_p   !< The cell locator object
    end subroutine set_volume

    !v Set the normal of specified face
    !
    !  Normalises the stored normal.
    module subroutine set_normal(loc_f, normal)
      type(face_locator), intent(in) :: loc_f            !< The face locator object
      real(ccs_real), dimension(:), intent(in) :: normal !< Array holding the face normal
    end subroutine set_normal

    !v Set face interpolation from cell and its local face id
    module subroutine set_face_interpolation(interpol_factor, loc_f)
      real(ccs_real), intent(in) :: interpol_factor  !< the interpolation factor to be used for loc_p
      type(face_locator), intent(inout) :: loc_f        !< the face locator object linking a cell-relative
    end subroutine set_face_interpolation

    !v Retrieves face interpolation from a face locator
    pure module subroutine get_face_interpolation(loc_f, interpol_factor)
      type(face_locator), intent(in) :: loc_f        !< the face locator object
      real(ccs_real), intent(out) :: interpol_factor  !< the interpolation factor to be used for loc_f
    end subroutine get_face_interpolation

    !> Query whether mesh was generated or read
    pure module subroutine get_mesh_generated(is_generated)
      logical, intent(out) :: is_generated !< The generated/read (true/false) status
    end subroutine

    !> Set whether a mesh was generated or read
    module subroutine set_mesh_generated(is_generated)
      logical, intent(in) :: is_generated
    end subroutine

    !> Get the numerical ID of a boundary from its name
    pure module subroutine get_bc_id(mesh, name, bc_id)
      type(ccs_mesh), intent(in) :: mesh     !< The mesh
      character(len=*), intent(in) :: name   !< The boundary name
      integer(ccs_int), intent(out) :: bc_id !< The boundary ID
    end subroutine get_bc_id

    !> Find a list of (boundary) faces based on their boundary name
    module subroutine find_face_entities(mesh, name, faces)
      type(ccs_mesh), intent(in) :: mesh                                  !< The mesh
      character(len=*), intent(in) :: name                                !< The boundary name
      type(face_locator), dimension(:), allocatable, intent(out) :: faces !< The list of boundary faces
    end subroutine find_face_entities
    
  end interface

end module meshing

