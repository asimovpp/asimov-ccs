submodule(meshing) meshing_accessors
#include "ccs_macros.inc"

  use utils, only: exit_print, str

  implicit none

contains

  !v Constructs a face locator object.
  !
  !  Creates the association between a face relative to a cell, i.e. to access the
  !  nth face of cell i.
  module subroutine set_face_location(mesh, index_p, cell_face_ctr, loc_f)
    type(ccs_mesh), target, intent(in) :: mesh      !< the mesh object being referred to.
    integer(ccs_int), intent(in) :: index_p         !< the index of the cell whose face is being accessed.
    integer(ccs_int), intent(in) :: cell_face_ctr   !< the cell-local index of the face.
    type(face_locator), intent(out) :: loc_f        !< the face locator object linking a cell-relative
    !< index with the mesh.
    loc_f%mesh => mesh
    loc_f%index_p = index_p
    loc_f%cell_face_ctr = cell_face_ctr
  end subroutine set_face_location

  !v Constructs a cell locator object.
  !
  !  Creates the association between a mesh and cell index, storing it in the
  !  returned cell locator object.
  module subroutine set_cell_location(mesh, index_p, loc_p)
    type(ccs_mesh), target, intent(in) :: mesh !< the mesh object being referred to.
    integer(ccs_int), intent(in) :: index_p    !< the cell index.
    type(cell_locator), intent(out) :: loc_p   !< the cell locator object linking a cell index with the mesh.

    loc_p%mesh => mesh
    loc_p%index_p = index_p

    ! XXX: Potentially expensive...
    if (index_p > mesh%topo%total_num_cells) then
      call error_abort("ERROR: trying to access cell I don't have access to." // str(index_p) // " " // str(mesh%topo%local_num_cells))
    end if
  end subroutine set_cell_location

  !v Constructs a neighbour locator object.
  !
  !  Creates the association between a neighbour cell F relative to cell P, i.e. to
  !  access the nth neighbour of cell i.
  module subroutine set_neighbour_location(loc_p, nb_counter, loc_nb)
    type(cell_locator), intent(in) :: loc_p        !< the cell locator object of the cell
    !< whose neighbour is being accessed.
    integer(ccs_int), intent(in) :: nb_counter     !< the cell-local index of the neighbour.
    type(neighbour_locator), intent(out) :: loc_nb !< the neighbour locator object linking a
    !< cell-relative index with the mesh.

    loc_nb%mesh => loc_p%mesh
    loc_nb%index_p = loc_p%index_p

    ! XXX: Safe, but would create a circular dependency...
    ! XXX: Potentially expensive...
    ! call count_neighbours(loc_p, nnb)
    ! if (nb_counter > nnb) then
    !   call error_abort("ERROR: cell has fewer neighbours than neighbour count requested.")
    ! else if (nb_counter < 1) then
    !   call error_abort("ERROR: cell neighbour counter must be >= 1.")
    ! else
    !   loc_nb%nb_counter = nb_counter
    ! end if

    loc_nb%nb_counter = nb_counter

    associate (mymesh => loc_nb%mesh, &
               i => loc_nb%index_p, &
               j => loc_nb%nb_counter)
      if (mymesh%topo%nb_indices(j, i) == i) then
        call error_abort("ERROR: attempted to set self as neighbour. Cell: " // str(i) // str(j))
      end if
    end associate
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

    loc_v%mesh => mesh
    loc_v%index_p = index_p
    loc_v%cell_vert_ctr = cell_vert_ctr
  end subroutine set_vert_location

  !> Set face index
  module subroutine set_face_index(index_p, cell_face_ctr, index_f, mesh)
    integer(ccs_int), intent(in) :: index_p
    integer(ccs_int), intent(in) :: cell_face_ctr
    integer(ccs_int), intent(in) :: index_f
    type(ccs_mesh), target, intent(inout) :: mesh

    mesh%topo%face_indices(cell_face_ctr, index_p) = index_f
  end subroutine set_face_index

  !> Returns the normal vector of a face
  module subroutine get_face_normal(loc_f, normal)
    type(face_locator), intent(in) :: loc_f                !< the face locator object.
    real(ccs_real), dimension(ndim), intent(out) :: normal !< an ndimensional array representing the face normal vector.

    associate (mesh => loc_f%mesh, &
               cell => loc_f%index_p, &
               face => loc_f%cell_face_ctr)
      normal(:) = mesh%geo%face_normals(:, face, cell)
    end associate
  end subroutine get_face_normal

  !> Returns the area of a face
  module subroutine get_face_area(loc_f, area)
    type(face_locator), intent(in) :: loc_f !< the face locator object.
    real(ccs_real), intent(out) :: area     !< the face area.

    associate (mesh => loc_f%mesh, &
               cell => loc_f%index_p, &
               face => loc_f%cell_face_ctr)
      area = mesh%geo%face_areas(face, cell)
    end associate
  end subroutine get_face_area

  !> Set the area of specified face
  module subroutine set_area(area, loc_f)
    real(ccs_real), intent(in) :: area      !< The face area
    type(face_locator), intent(in) :: loc_f !< The face locator object

    associate(mesh => loc_f%mesh, &
         cell => loc_f%index_p, &
         face => loc_f%cell_face_ctr)
      mesh%geo%face_areas(face, cell) = area
    end associate
  end subroutine set_area

  !> Returns the centre of a cell
  module subroutine get_cell_centre(loc_p, x)
    type(cell_locator), intent(in) :: loc_p           !< the cell locator object.
    real(ccs_real), dimension(:), intent(out) :: x !< an ndimensional array representing the cell centre.

    integer :: dim
    
    associate (mesh => loc_p%mesh, &
               cell => loc_p%index_p)
      do dim = 1, min(size(x), ndim)
         x(dim) = mesh%geo%x_p(dim, cell)
      end do
    end associate
  end subroutine get_cell_centre

  !> Returns the centre of a neighbour cell
  module subroutine get_neighbour_centre(loc_nb, x)
    type(neighbour_locator), intent(in) :: loc_nb     !< the neighbour locator object.
    real(ccs_real), dimension(ndim), intent(out) :: x !< an ndimensional array representing the neighbour cell centre.

    type(cell_locator) :: cell_loc_nb

    call get_neighbour_cell_locator(loc_nb, cell_loc_nb)
    call get_cell_centre(cell_loc_nb, x)
  end subroutine get_neighbour_centre

  !> Returns the centre of a face
  module subroutine get_face_centre(loc_f, x)
    type(face_locator), intent(in) :: loc_f           !< the face locator object.
    real(ccs_real), dimension(ndim), intent(out) :: x !< an ndimensional array representing the face centre.

    associate (mesh => loc_f%mesh, &
               cell => loc_f%index_p, &
               face => loc_f%cell_face_ctr)
      x(:) = mesh%geo%x_f(:, face, cell)
    end associate
  end subroutine get_face_centre

  !> Returns the centre of a vertex
  module subroutine get_vert_centre(loc_v, x)
    type(vert_locator), intent(in) :: loc_v           !< the vertex locator object.
    real(ccs_real), dimension(:), intent(out) :: x !< an ndimensional array representing the vertex centre.

    integer :: dim

    associate(mesh => loc_v%mesh, &
      cell => loc_v%index_p, &
      vert => loc_v%cell_vert_ctr)
      do dim = 1, min(size(x), ndim)
         x(dim) = mesh%geo%vert_coords(dim, vert, cell)
      end do
    end associate
  end subroutine get_vert_centre

  !> Returns the volume of a cell
  module subroutine get_cell_volume(loc_p, V)
    type(cell_locator), intent(in) :: loc_p !< the cell locator object.
    real(ccs_real), intent(out) :: V        !< the cell volume.

    associate (mesh => loc_p%mesh, &
               cell => loc_p%index_p)
      V = mesh%geo%volumes(cell)
    end associate
  end subroutine get_cell_volume

  !> Returns the volume of a neighbour cell
  module subroutine get_neighbour_volume(loc_nb, V)
    type(neighbour_locator), intent(in) :: loc_nb !< the neighbour locator object.
    real(ccs_real), intent(out) :: V              !< the neighbour cell volume.

    type(cell_locator) :: cell_loc_nb

    call get_neighbour_cell_locator(loc_nb, cell_loc_nb)
    call get_cell_volume(cell_loc_nb, V)
  end subroutine get_neighbour_volume

  !> Returns the global index of a cell
  module subroutine get_cell_global_index(loc_p, global_index_p)
    type(cell_locator), intent(in) :: loc_p         !< the cell locator object.
    integer(ccs_int), intent(out) :: global_index_p !< the global index of the cell.

    associate (mesh => loc_p%mesh)
      if (mesh%topo%local_num_cells > 0) then ! XXX: Potentially expensive...
        associate (cell => loc_p%index_p)
          global_index_p = mesh%topo%global_indices(cell)
        end associate
      else
        global_index_p = -1 ! XXX: What should we do in case of too many processors for a given mesh?
      end if
    end associate
  end subroutine get_cell_global_index

  !> Returns the global index of a neighbouring cell
  module subroutine get_neighbour_global_index(loc_nb, global_index_nb)
    type(neighbour_locator), intent(in) :: loc_nb    !< the neighbour locator object.
    integer(ccs_int), intent(out) :: global_index_nb !< the global index of the neighbour cell.

    type(cell_locator) :: cell_loc_nb
    call get_neighbour_cell_locator(loc_nb, cell_loc_nb)
    call get_global_index(cell_loc_nb, global_index_nb)
  end subroutine get_neighbour_global_index

  !> Returns the neighbour count of a cell (including boundary neighbours)
  module subroutine cell_count_neighbours(loc_p, nnb)
    type(cell_locator), intent(in) :: loc_p !< the cell locator object.
    integer(ccs_int), intent(out) :: nnb    !< the neighbour count of the cell.

    associate (mesh => loc_p%mesh, &
               cell => loc_p%index_p)
      nnb = mesh%topo%num_nb(cell)
    end associate
  end subroutine cell_count_neighbours

  !> Returns the boundary status of a neighbouring cell
  module subroutine get_neighbour_boundary_status(loc_nb, is_boundary)
    type(neighbour_locator), intent(in) :: loc_nb !< the neighbour locator object.
    logical, intent(out) :: is_boundary           !< the boundary status of the neighbour.

    integer :: index_nb

    call get_neighbour_local_index(loc_nb, index_nb)

    if (index_nb > 0) then
      is_boundary = .false.
    else if (index_nb < 0) then
      is_boundary = .true.
    else
      call error_abort("ERROR: neighbour index (0) is invalid.")
    end if
  end subroutine get_neighbour_boundary_status

  !> Returns the boundary status of a face
  module subroutine get_face_boundary_status(loc_f, is_boundary)
    type(face_locator), intent(in) :: loc_f !< the face locator object.
    logical, intent(out) :: is_boundary     !< the boundary status of the neighbour.

    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb

    associate (mesh => loc_f%mesh, &
               i => loc_f%index_p, &
               j => loc_f%cell_face_ctr)
      call set_cell_location(mesh, i, loc_p)
      call set_neighbour_location(loc_p, j, loc_nb)
    end associate
    call get_neighbour_boundary_status(loc_nb, is_boundary)
  end subroutine get_face_boundary_status

  !v Returns the local distribution status of a neighbouring cell
  !
  !  Given a distributed mesh, a processor needs both the cells within its partition
  !  and cells from the surrounding halo - this subroutine get_indicates whether a
  !  cell's neighbour is within the local partition or the halo.
  module subroutine get_local_status(loc_nb, is_local)
    type(neighbour_locator), intent(in) :: loc_nb !< the neighbour locator object.
    logical, intent(out) :: is_local !< the local status of the neighbour.

    integer :: index_nb

    call get_neighbour_local_index(loc_nb, index_nb)
    associate (mesh => loc_nb%mesh)
      if ((index_nb > 0) .and. (index_nb <= mesh%topo%local_num_cells)) then
        is_local = .true.
      else
        is_local = .false.
      end if
    end associate
  end subroutine get_local_status

  !v Returns the local index of a cell
  !
  !  Generally the local index of a cell is should be the same as its location within
  !  the local cell vector - this particular subroutine get_is therefore expected of
  !  limited use and is mostly present for uniformity.
  module subroutine get_cell_local_index(loc_p, index_p)
    type(cell_locator), intent(in) :: loc_p  !< the cell locator object.
    integer(ccs_int), intent(out) :: index_p !< the local index of the cell.

    index_p = loc_p%index_p
  end subroutine get_cell_local_index

  !> Returns the local index of a neighbouring cell
  module subroutine get_neighbour_local_index(loc_nb, index_nb)
    type(neighbour_locator), intent(in) :: loc_nb !< the neighbour locator object.
    integer(ccs_int), intent(out) :: index_nb     !< the local index of the neighbour cell.

    associate (mesh => loc_nb%mesh, &
               i => loc_nb%index_p, &
               j => loc_nb%nb_counter)
      index_nb = mesh%topo%nb_indices(j, i)
    end associate
  end subroutine get_neighbour_local_index

  !> Returns the local index of a face
  module subroutine get_face_local_index(loc_f, index_f)
    type(face_locator), intent(in) :: loc_f  !< the face locator object.
    integer(ccs_int), intent(out) :: index_f !< the local index of the face.

    associate (mesh => loc_f%mesh, &
               i => loc_f%index_p, &
               j => loc_f%cell_face_ctr)
      index_f = mesh%topo%face_indices(j, i)
    end associate
  end subroutine get_face_local_index

  subroutine get_neighbour_cell_locator(loc_nb, loc_p)
    type(neighbour_locator), intent(in) :: loc_nb
    type(cell_locator), intent(out) :: loc_p

    integer(ccs_int) :: index_nb

    call get_local_index(loc_nb, index_nb)
    call set_cell_location(loc_nb%mesh, index_nb, loc_p)
  end subroutine get_neighbour_cell_locator

  !> Set the cell centre of specified cell
  module subroutine set_cell_centre(loc_p, x_p)
    type(cell_locator), intent(in) :: loc_p         !< The cell locator object.
    real(ccs_real), dimension(:), intent(in) :: x_p !< The cell centre array.

    integer :: dim

    associate(mesh => loc_p%mesh, &
         i => loc_p%index_p)
      do dim = 1, min(size(x_p), ndim)
         mesh%geo%x_p(dim, i) = x_p(dim)
      end do
    end associate
  end subroutine set_cell_centre

  !> Set the face centre of specified face
  module subroutine set_face_centre(loc_f, x_f)
    type(face_locator), intent(in) :: loc_f         !< The face locator object.
    real(ccs_real), dimension(:), intent(in) :: x_f !< The face centre array.

    integer :: dim

    associate(mesh => loc_f%mesh, &
         i => loc_f%index_p, &
         j => loc_f%cell_face_ctr)
      do dim = 1, min(size(x_f), ndim)
         mesh%geo%x_f(dim, j, i) = x_f(dim)
      end do
    end associate
  end subroutine set_face_centre

  !> Set the centre of specified vertex
  module subroutine set_vert_centre(loc_v, x_v)
    type(vert_locator), intent(in) :: loc_v         !< The vertex locator object.
    real(ccs_real), dimension(:), intent(in) :: x_v !< The vertex centre array.

    integer :: dim

    associate(mesh => loc_v%mesh, &
         i => loc_v%index_p, &
         j => loc_v%cell_vert_ctr)
      do dim = 1, min(size(x_v), ndim)
         mesh%geo%vert_coords(dim, j, i) = x_v(dim)
      end do
    end associate
  end subroutine set_vert_centre

  !v Set the normal of specified face
  !
  !  Normalises the stored normal.
  module subroutine set_normal(loc_f, normal)
    type(face_locator), intent(in) :: loc_f            !< The face locator object
    real(ccs_real), dimension(:), intent(in) :: normal !< Array holding the face normal

    integer :: dim
    real(ccs_real) :: invmag

    invmag = 1.0_ccs_real / sqrt(sum(normal**2))
    associate(mesh => loc_f%mesh, &
         cell => loc_f%index_p, &
         face => loc_f%cell_face_ctr)
      do dim = 1, min(size(normal), ndim)
         mesh%geo%face_normals(dim, face, cell) = normal(dim) * invmag
      end do
    end associate
  end subroutine set_normal
  
end submodule meshing_accessors
