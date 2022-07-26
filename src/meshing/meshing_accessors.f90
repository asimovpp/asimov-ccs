submodule(meshing) meshing_accessors
#include "ccs_macros.inc"

  use utils, only: exit_print, str

  implicit none

contains

  !>  Constructs a face locator object.
  !
  !> @description Creates the association between a face relative to a cell, i.e. to access the
  !!              nth face of cell i.
  !
  !> @param[in]  mesh         mesh
  !> @param[in]  ccs_int     index_p
  !> @param[in]  ccs_int     cell_face_ctr
  !> @param[out] face_locator loc_f
  !!
  module subroutine set_face_location(mesh, index_p, cell_face_ctr, loc_f)
    type(ccs_mesh), target, intent(in) :: mesh      !< the mesh object being referred to.
    integer(ccs_int), intent(in) :: index_p         !< the index of the cell whose face is being accessed.
    integer(ccs_int), intent(in) :: cell_face_ctr   !< the cell-local index of the face.
    type(face_locator), intent(out) :: loc_f        !< the face locator object linking a cell-relative
                                                    !!  index with the mesh.
    loc_f%mesh => mesh
    loc_f%index_p = index_p
    loc_f%cell_face_ctr = cell_face_ctr
  end subroutine set_face_location

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

    loc_p%mesh => mesh
    loc_p%index_p = index_p

    ! XXX: Potentially expensive...
    if (index_p > mesh%ntotal) then
      call error_abort("ERROR: trying to access cell I don't have access to!" // str(index_p) // str(mesh%nlocal))
    end if
  end subroutine set_cell_location

  !>  Constructs a neighbour locator object.
  !
  !> @description Creates the association between a neighbour cell F relative to cell P, i.e. to
  !!              access the nth neighbour of cell i.
  !
  !> @param[in]  cell_locator      loc_p      - the cell locator object of the cell whose
  !!                                                    neighbour is being accessed.
  !> @param[in]  ccs_int           nb_counter - the cell-local index of the neighbour.
  !> @param[out] neighbour_locator neighbour_location - the neighbour locator object linking a
  !!                                                    cell-relative index with the mesh.
  module subroutine set_neighbour_location(loc_p, nb_counter, loc_nb)
    type(cell_locator), intent(in) :: loc_p
    integer(ccs_int), intent(in) :: nb_counter
    type(neighbour_locator), intent(out) :: loc_nb

    loc_nb%mesh => loc_p%mesh
    loc_nb%index_p = loc_p%index_p

    !! XXX: Safe, but would create a circular dependency...
    !! ! XXX: Potentially expensive...
    !! call count_neighbours(loc_p, nnb)
    !! if (nb_counter > nnb) then
    !!   call error_abort("ERROR: cell has fewer neighbours than neighbour count requested!")
    !! else if (nb_counter < 1) then
    !!   call error_abort("ERROR: cell neighbour counter must be >= 1!")
    !! else
    !!   loc_nb%nb_counter = nb_counter
    !! end if

    loc_nb%nb_counter = nb_counter

    associate (mymesh => loc_nb%mesh, &
               i => loc_nb%index_p, &
               j => loc_nb%nb_counter)
      if (mymesh%neighbour_indices(j, i) == i) then
        call error_abort("ERROR: trying to set self as neighbour! Cell: " // str(i) // str(j))
      end if
    end associate
  end subroutine set_neighbour_location

  !>  Set face index
  module subroutine set_face_index(index_p, cell_face_ctr, index_f, mesh)
    integer(ccs_int), intent(in) :: index_p
    integer(ccs_int), intent(in) :: cell_face_ctr
    integer(ccs_int), intent(in) :: index_f
    type(ccs_mesh), target, intent(inout) :: mesh

    mesh%face_indices(cell_face_ctr, index_p) = index_f
  end subroutine set_face_index

  !>  Returns the normal vector of a face
  !
  !> @param[in]  face_locator    loc_f - the face locator object.
  !> @param[out] real(ccs_real) normal(ndim)  - an ndimensional array representing the face normal
  !!                                             vector.
  module subroutine get_face_normal(loc_f, normal)
    type(face_locator), intent(in) :: loc_f
    real(ccs_real), dimension(ndim), intent(out) :: normal

    associate (mesh => loc_f%mesh, &
               cell => loc_f%index_p, &
               face => loc_f%cell_face_ctr)
      normal(:) = mesh%face_normals(:, face, cell)
    end associate
  end subroutine get_face_normal

  !>  Returns the area of a face
  !
  !> @param[in]  face_locator    loc_f - the face locator object.
  !> @param[out] real(ccs_real) area          - the face area.
  module subroutine get_face_area(loc_f, area)
    type(face_locator), intent(in) :: loc_f
    real(ccs_real), intent(out) :: area

    associate (mesh => loc_f%mesh, &
               cell => loc_f%index_p, &
               face => loc_f%cell_face_ctr)
      area = mesh%face_areas(face, cell)
    end associate
  end subroutine get_face_area

  !>  Returns the centre of a cell
  !
  !> @param[in]  cell_locator     loc_p - the cell locator object.
  !> @param[out] real(ccs_real)  x(ndim)       - an ndimensional array representing the cell centre.
  module subroutine get_cell_centre(loc_p, x)
    type(cell_locator), intent(in) :: loc_p
    real(ccs_real), dimension(ndim), intent(out) :: x

    associate (mesh => loc_p%mesh, &
               cell => loc_p%index_p)
      x(:) = mesh%x_p(:, cell)
    end associate
  end subroutine get_cell_centre

  !>  Returns the centre of a neighbour cell
  !
  !> @param[in]  neighbour_locator loc_nb - the neighbour locator object.
  !> @param[out] real(ccs_real)   x(ndim)            - an ndimensional array representing the
  !!                                                    neighbour cell centre.
  module subroutine get_neighbour_centre(loc_nb, x)
    type(neighbour_locator), intent(in) :: loc_nb
    real(ccs_real), dimension(ndim), intent(out) :: x

    type(cell_locator) :: cell_loc_nb

    call get_neighbour_cell_locator(loc_nb, cell_loc_nb)
    call get_cell_centre(cell_loc_nb, x)
  end subroutine get_neighbour_centre

  !>  Returns the centre of a face
  !
  !> @param[in]  face_locator     loc_f - the face locator object.
  !> @param[out] real(ccs_real)  x(ndim)       - an ndimensional array representing the face centre.
  module subroutine get_face_centre(loc_f, x)
    type(face_locator), intent(in) :: loc_f
    real(ccs_real), dimension(ndim), intent(out) :: x

    associate (mesh => loc_f%mesh, &
               cell => loc_f%index_p, &
               face => loc_f%cell_face_ctr)
      x(:) = mesh%x_f(:, face, cell)
    end associate
  end subroutine get_face_centre

  !>  Returns the volume of a cell
  !
  !> @param[in] cell_locator     loc_p - the cell locator object.
  !> @param[out] real(ccs_real) V             - the cell volume.
  module subroutine get_cell_volume(loc_p, V)
    type(cell_locator), intent(in) :: loc_p
    real(ccs_real), intent(out) :: V

    associate (mesh => loc_p%mesh, &
               cell => loc_p%index_p)
      V = mesh%volumes(cell)
    end associate
  end subroutine get_cell_volume

  !>  Returns the volume of a neighbour cell
  !
  !> @param[in] neighbour_locator loc_nb - the neighbour locator object.
  !> @param[out] real(ccs_real)  V                  - the neighbour cell volume.
  module subroutine get_neighbour_volume(loc_nb, V)
    type(neighbour_locator), intent(in) :: loc_nb
    real(ccs_real), intent(out) :: V

    type(cell_locator) :: cell_loc_nb

    call get_neighbour_cell_locator(loc_nb, cell_loc_nb)
    call get_cell_volume(cell_loc_nb, V)
  end subroutine get_neighbour_volume

  !>  Returns the global index of a cell
  !
  !> @param[in]  cell_locator      loc_p - the cell locator object.
  !> @param[out] integer(ccs_int) global_index_p          - the global index of the cell.
  module subroutine get_cell_global_index(loc_p, global_index_p)
    type(cell_locator), intent(in) :: loc_p
    integer(ccs_int), intent(out) :: global_index_p

    associate (mesh => loc_p%mesh)
      if (mesh%nlocal > 0) then ! XXX: Potentially expensive...
        associate (cell => loc_p%index_p)
          global_index_p = mesh%global_indices(cell)
        end associate
      else
        global_index_p = -1 ! XXX: What should we do in case of too many processors for a given mesh?
      end if
    end associate
  end subroutine get_cell_global_index

  !>  Returns the global index of a neighbouring cell
  !
  !> @param[in]  neighbour_locator loc_nb           - the neighbour locator object.
  !> @param[out] integer(ccs_int) global_index_nb   - the global index of the neighbour cell.
  module subroutine get_neighbour_global_index(loc_nb, global_index_nb)
    type(neighbour_locator), intent(in) :: loc_nb
    integer(ccs_int), intent(out) :: global_index_nb

    type(cell_locator) :: cell_loc_nb
    call get_neighbour_cell_locator(loc_nb, cell_loc_nb)
    call get_global_index(cell_loc_nb, global_index_nb)
  end subroutine get_neighbour_global_index

  !>  Returns the neighbour count of a cell (including boundary neighbours)
  !
  !> @param[in]  cell_locator      loc_p - the cell locator object.
  !> @param[out] integer(ccs_int) nnb           - the neighbour count of the cell.
  module subroutine cell_count_neighbours(loc_p, nnb)
    type(cell_locator), intent(in) :: loc_p
    integer(ccs_int), intent(out) :: nnb

    associate (mesh => loc_p%mesh, &
               cell => loc_p%index_p)
      nnb = mesh%nnb(cell)
    end associate
  end subroutine cell_count_neighbours

  !>  Returns the boundary status of a neighbouring cell
  !
  !> @param[in]  neighbour_locator loc_nb - the neighbour locator object.
  !> @param[out] logical           is_boundary        - the boundary status of the neighbour.
  module subroutine get_neighbour_boundary_status(loc_nb, is_boundary)
    type(neighbour_locator), intent(in) :: loc_nb
    logical, intent(out) :: is_boundary

    integer :: index_nb

    call get_neighbour_local_index(loc_nb, index_nb)

    if (index_nb > 0) then
      is_boundary = .false.
    else if (index_nb < 0) then
      is_boundary = .true.
    else
      call error_abort("ERROR: neighbour index (0) is invalid")
    end if
  end subroutine get_neighbour_boundary_status

  !>  Returns the boundary status of a face
  !
  !> @param[in]  face_locator loc_f - the face locator object.
  !> @param[out] logical      is_boundary   - the boundary status of the neighbour.
  module subroutine get_face_boundary_status(loc_f, is_boundary)
    type(face_locator), intent(in) :: loc_f
    logical, intent(out) :: is_boundary

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

    integer :: index_nb

    call get_neighbour_local_index(loc_nb, index_nb)
    associate (mesh => loc_nb%mesh)
      if ((index_nb > 0) .and. (index_nb <= mesh%nlocal)) then
        is_local = .true.
      else
        is_local = .false.
      end if
    end associate
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

    index_p = loc_p%index_p
  end subroutine get_cell_local_index

  !>  Returns the local index of a neighbouring cell
  !
  !> @param[in]  neighbour_locator loc_nb - the neighbour locator object.
  !> @param[out] integer(ccs_int) index_nb              - the local index of the neighbour cell.
  module subroutine get_neighbour_local_index(loc_nb, index_nb)
    type(neighbour_locator), intent(in) :: loc_nb
    integer(ccs_int), intent(out) :: index_nb

    associate (mesh => loc_nb%mesh, &
               i => loc_nb%index_p, &
               j => loc_nb%nb_counter)
      index_nb = mesh%neighbour_indices(j, i)
    end associate
  end subroutine get_neighbour_local_index

  !>  Returns the local index of a face
  !
  !> @param[in]  face_locator     loc_f    - the face locator object.
  !> @param[out] integer(ccs_int) index_f  - the local index of the face.
  module subroutine get_face_local_index(loc_f, index_f)
    type(face_locator), intent(in) :: loc_f
    integer(ccs_int), intent(out) :: index_f

    associate (mesh => loc_f%mesh, &
               i => loc_f%index_p, &
               j => loc_f%cell_face_ctr)
      index_f = mesh%face_indices(j, i)
    end associate
  end subroutine get_face_local_index

  subroutine get_neighbour_cell_locator(loc_nb, loc_p)
    type(neighbour_locator), intent(in) :: loc_nb
    type(cell_locator), intent(out) :: loc_p

    integer(ccs_int) :: index_nb

    call get_local_index(loc_nb, index_nb)
    call set_cell_location(loc_nb%mesh, index_nb, loc_p)
  end subroutine get_neighbour_cell_locator

end submodule meshing_accessors
