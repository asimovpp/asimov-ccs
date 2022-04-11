submodule (meshing) meshing_accessors

  implicit none
  
contains
  
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
    type(ccs_mesh), target, intent(in) :: geometry
    integer(ccs_int), intent(in) :: cell_idx
    integer(ccs_int), intent(in) :: cell_face_ctr
    type(face_locator), intent(out) :: face_location

    face_location%mesh => geometry
    face_location%cell_idx = cell_idx
    face_location%cell_face_ctr = cell_face_ctr
  end subroutine set_face_location

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
    type(ccs_mesh), target, intent(in) :: geometry
    integer(ccs_int), intent(in) :: cell_idx
    type(cell_locator), intent(out) :: cell_location

    ! XXX: Potentially expensive...
    if (cell_idx > geometry%ntotal) then
      print *, "ERROR: trying to access cell I don't have access to!", cell_idx, geometry%nlocal
      stop
    else
      cell_location%mesh => geometry
      cell_location%cell_idx = cell_idx
    end if
  end subroutine set_cell_location
  
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

  !> @brief Set face index
  module subroutine set_face_index(cell_idx, cell_face_ctr, face_idx, geometry)
    integer(ccs_int), intent(in) :: cell_idx
    integer(ccs_int), intent(in) :: cell_face_ctr
    integer(ccs_int), intent(in) :: face_idx
    type(ccs_mesh), target, intent(inout) :: geometry

    geometry%faceidx(cell_face_ctr, cell_idx) = face_idx
  end subroutine set_face_index

  !> @brief Returns the normal vector of a face
  !
  !> @param[in]  face_locator    face_location - the face locator object.
  !> @param[out] real(ccs_real) normal(ndim)  - an ndimensional array representing the face normal
  !!                                             vector.
  module subroutine get_face_normal(face_location, normal)
    type(face_locator), intent(in) :: face_location
    real(ccs_real), dimension(ndim), intent(out) :: normal

    associate(mesh => face_location%mesh, &
         cell => face_location%cell_idx, &
         face => face_location%cell_face_ctr)
      normal(:) = mesh%nf(:, face, cell)
    end associate
  end subroutine get_face_normal

  !> @brief Returns the area of a face
  !
  !> @param[in]  face_locator    face_location - the face locator object.
  !> @param[out] real(ccs_real) area          - the face area.
  module subroutine get_face_area(face_location, area)
    type(face_locator), intent(in) :: face_location
    real(ccs_real), intent(out) :: area

    associate(mesh => face_location%mesh, &
         cell =>face_location%cell_idx, &
         face =>face_location%cell_face_ctr)
      area = mesh%Af(face, cell)
    end associate
  end subroutine get_face_area

  !> @brief Returns the centre of a cell
  !
  !> @param[in]  cell_locator     cell_location - the cell locator object.
  !> @param[out] real(ccs_real)  x(ndim)       - an ndimensional array representing the cell centre.
  module subroutine get_cell_centre(cell_location, x)
    type(cell_locator), intent(in) :: cell_location
    real(ccs_real), dimension(ndim), intent(out) :: x

    associate(mesh => cell_location%mesh, &
         cell => cell_location%cell_idx)
      x(:) = mesh%xc(:, cell)
    end associate
  end subroutine get_cell_centre

  !> @brief Returns the centre of a neighbour cell
  !
  !> @param[in]  neighbour_locator neighbour_location - the neighbour locator object.
  !> @param[out] real(ccs_real)   x(ndim)            - an ndimensional array representing the
  !!                                                    neighbour cell centre.
  module subroutine get_neighbour_centre(neighbour_location, x)
    type(neighbour_locator), intent(in) :: neighbour_location
    real(ccs_real), dimension(ndim), intent(out) :: x

    type(cell_locator) :: nb_cell_location
    
    call get_neighbour_cell_locator(neighbour_location, nb_cell_location)
    call get_cell_centre(nb_cell_location, x)
  end subroutine get_neighbour_centre
  
  !> @brief Returns the centre of a face
  !
  !> @param[in]  face_locator     face_location - the face locator object.
  !> @param[out] real(ccs_real)  x(ndim)       - an ndimensional array representing the face centre.
  module subroutine get_face_centre(face_location, x)
    type(face_locator), intent(in) :: face_location
    real(ccs_real), dimension(ndim), intent(out) :: x

    associate(mesh => face_location%mesh, &
         cell => face_location%cell_idx, &
         face => face_location%cell_face_ctr)
      x(:) = mesh%xf(:, face, cell)
    end associate
  end subroutine get_face_centre

  !> @brief Returns the volume of a cell
  !
  !> @param[in] cell_locator     cell_location - the cell locator object.
  !> @param[out] real(ccs_real) V             - the cell volume.
  module subroutine get_cell_volume(cell_location, V)
    type(cell_locator), intent(in) :: cell_location
    real(ccs_real), intent(out) :: V

    associate(mesh => cell_location%mesh, &
         cell => cell_location%cell_idx)
      V = mesh%vol(cell)
    end associate
  end subroutine get_cell_volume

  !> @brief Returns the volume of a neighbour cell
  !
  !> @param[in] neighbour_locator neighbour_location - the neighbour locator object.
  !> @param[out] real(ccs_real)  V                  - the neighbour cell volume.
  module subroutine get_neighbour_volume(neighbour_location, V)
    type(neighbour_locator), intent(in) :: neighbour_location
    real(ccs_real), intent(out) :: V

    type(cell_locator) :: nb_cell_location

    call get_neighbour_cell_locator(neighbour_location, nb_cell_location)
    call get_cell_volume(nb_cell_location, V)
  end subroutine get_neighbour_volume

  !> @brief Returns the global index of a cell
  !
  !> @param[in]  cell_locator      cell_location - the cell locator object.
  !> @param[out] integer(ccs_int) idxg          - the global index of the cell.
  module subroutine get_cell_global_index(cell_location, idxg)
    type(cell_locator), intent(in) :: cell_location
    integer(ccs_int), intent(out) :: idxg

    associate(mesh => cell_location%mesh)
      if (mesh%nlocal > 0) then ! XXX: Potentially expensive...
        associate(cell => cell_location%cell_idx)
          idxg = mesh%idx_global(cell)
        end associate
      else
        idxg = -1 ! XXX: What should we do in case of too many processors for a given mesh?
      end if
    end associate
  end subroutine get_cell_global_index

  !> @brief Returns the global index of a neighbouring cell
  !
  !> @param[in]  neighbour_locator neighbour_location - the neighbour locator object.
  !> @param[out] integer(ccs_int) nbidxg             - the global index of the neighbour cell.
  module subroutine get_neighbour_global_index(neighbour_location, nbidxg)
    type(neighbour_locator), intent(in) :: neighbour_location
    integer(ccs_int), intent(out) :: nbidxg

    type(cell_locator) :: nb_cell_location
    call get_neighbour_cell_locator(neighbour_location, nb_cell_location)
    call get_global_index(nb_cell_location, nbidxg)
  end subroutine get_neighbour_global_index

  !> @brief Returns the neighbour count of a cell (including boundary neighbours)
  !
  !> @param[in]  cell_locator      cell_location - the cell locator object.
  !> @param[out] integer(ccs_int) nnb           - the neighbour count of the cell.
  module subroutine cell_count_neighbours(cell_location, nnb)
    type(cell_locator), intent(in) :: cell_location
    integer(ccs_int), intent(out) :: nnb

    associate(mesh => cell_location%mesh, &
         cell => cell_location%cell_idx)
      nnb = mesh%nnb(cell)
    end associate
  end subroutine cell_count_neighbours

  !> @brief Returns the boundary status of a neighbouring cell
  !
  !> @param[in]  neighbour_locator neighbour_location - the neighbour locator object.
  !> @param[out] logical           is_boundary        - the boundary status of the neighbour.
  module subroutine get_neighbour_boundary_status(neighbour_location, is_boundary)
    type(neighbour_locator), intent(in) :: neighbour_location
    logical, intent(out) :: is_boundary

    integer :: nbidx

    call get_neighbour_local_index(neighbour_location, nbidx)

    if (nbidx > 0) then
      is_boundary = .false.
    else if (nbidx < 0) then
      is_boundary = .true.
    else
      print *, "ERROR: neighbour index (0) is invalid"
      stop
    end if
  end subroutine get_neighbour_boundary_status

  !> @brief Returns the boundary status of a face
  !
  !> @param[in]  face_locator face_location - the face locator object.
  !> @param[out] logical      is_boundary   - the boundary status of the neighbour.
  module subroutine get_face_boundary_status(face_location, is_boundary)
    type(face_locator), intent(in) :: face_location
    logical, intent(out) :: is_boundary

    type(cell_locator) :: cell_location
    type(neighbour_locator) :: neighbour_location

    associate(mesh => face_location%mesh, &
         i => face_location%cell_idx, &
         j => face_location%cell_face_ctr)
      call set_cell_location(mesh, i, cell_location)
      call set_neighbour_location(cell_location, j, neighbour_location)
    end associate
    call get_neighbour_boundary_status(neighbour_location, is_boundary)
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
    
    integer :: nbidx

    call get_neighbour_local_index(neighbour_location, nbidx)
    associate(mesh => neighbour_location%mesh)
      if ((nbidx > 0) .and. (nbidx <= mesh%nlocal)) then
        is_local = .true.
      else
        is_local = .false.
      end if
    end associate
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
    
    idx = cell_location%cell_idx
  end subroutine get_cell_local_index

  !> @brief Returns the local index of a neighbouring cell
  !
  !> @param[in]  neighbour_locator neighbour_location - the neighbour locator object.
  !> @param[out] integer(ccs_int) nbidx              - the local index of the neighbour cell.
  module subroutine get_neighbour_local_index(neighbour_location, nbidx)
    type(neighbour_locator), intent(in) :: neighbour_location
    integer(ccs_int), intent(out) :: nbidx

    associate(mesh => neighbour_location%mesh, &
         i => neighbour_location%cell_idx, &
         j => neighbour_location%cell_neighbour_ctr)
      nbidx = mesh%nbidx(j, i)
    end associate
  end subroutine get_neighbour_local_index

  !> @brief Returns the local index of a face
  !
  !> @param[in]  face_locator      face_location - the face locator object.
  !> @param[out] integer(ccs_int) idx           - the local index of the face.
  module subroutine get_face_local_index(face_location, idx)
    type(face_locator), intent(in) :: face_location
    integer(ccs_int), intent(out) :: idx

    associate(mesh => face_location%mesh, &
      i => face_location%cell_idx, &
      j => face_location%cell_face_ctr)
      idx = mesh%faceidx(j,i)
    end associate
  end subroutine get_face_local_index

  subroutine get_neighbour_cell_locator(neighbour_location, cell_location)
    type(neighbour_locator), intent(in) :: neighbour_location
    type(cell_locator), intent(out) :: cell_location

    integer(ccs_int) :: nbidx
    
    call get_local_index(neighbour_location, nbidx)
    call set_cell_location(neighbour_location%mesh, nbidx, cell_location)
  end subroutine get_neighbour_cell_locator

end submodule meshing_accessors
