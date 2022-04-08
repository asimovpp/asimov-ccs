submodule (meshing) meshing_accessors

  implicit none
  
contains
  
  !> @brief Constructs a face locator object.
  !
  !> @description Creates the association between a face relative to a cell, i.e. to access the
  !!              nth face of cell i.
  !
  !> @param[in]  mesh         geometry      - the mesh object being referred to.
  !> @param[in]  accs_int     cell_idx      - the index of the cell whose face is being accessed.
  !> @param[in]  accs_int     cell_face_ctr - the cell-local index of the face.
  !> @param[out] face_locator face_location - the face locator object linking a cell-relative
  !!                                          index with the mesh.
  module subroutine set_face_location(geometry, cell_idx, cell_face_ctr, face_location)
    type(mesh), target, intent(in) :: geometry
    integer(accs_int), intent(in) :: cell_idx
    integer(accs_int), intent(in) :: cell_face_ctr
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
  !> @param[in]  accs_int     cell_idx      - the cell index. 
  !> @param[out] cell_locator cell_location - the cell locator object linking a cell index with
  !!                                          the mesh.
  module subroutine set_cell_location(geometry, cell_idx, cell_location)
    type(mesh), target, intent(in) :: geometry
    integer(accs_int), intent(in) :: cell_idx
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
  !> @param[in]  accs_int          cell_neighbour_ctr - the cell-local index of the neighbour.
  !> @param[out] neighbour_locator neighbour_location - the neighbour locator object linking a
  !!                                                    cell-relative index with the mesh.
  module subroutine set_neighbour_location(cell_location, cell_neighbour_ctr, neighbour_location)
    type(cell_locator), intent(in) :: cell_location
    integer(accs_int), intent(in) :: cell_neighbour_ctr
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

  module procedure set_face_index
    geometry%faceidx(cell_face_ctr, cell_idx) = face_idx

  end procedure set_face_index


  module procedure get_face_normal
    associate(mesh => face_location%mesh, &
         cell => face_location%cell_idx, &
         face => face_location%cell_face_ctr)
      normal(:) = mesh%nf(:, face, cell)
    end associate
  end procedure get_face_normal

  module procedure get_face_area
    associate(mesh => face_location%mesh, &
         cell =>face_location%cell_idx, &
         face =>face_location%cell_face_ctr)
      area = mesh%Af(face, cell)
    end associate
  end procedure get_face_area

  module procedure get_cell_centre
    associate(mesh => cell_location%mesh, &
         cell => cell_location%cell_idx)
      x(:) = mesh%xc(:, cell)
    end associate
  end procedure get_cell_centre

  module procedure get_neighbour_centre
    type(cell_locator) :: nb_cell_location
    
    call get_neighbour_cell_locator(neighbour_location, nb_cell_location)
    call get_cell_centre(nb_cell_location, x)
  end procedure get_neighbour_centre
  
  module procedure get_face_centre
    associate(mesh => face_location%mesh, &
         cell => face_location%cell_idx, &
         face => face_location%cell_face_ctr)
      x(:) = mesh%xf(:, face, cell)
    end associate
  end procedure get_face_centre

  module procedure get_cell_volume
    associate(mesh => cell_location%mesh, &
         cell => cell_location%cell_idx)
      V = mesh%vol(cell)
    end associate
  end procedure get_cell_volume

  module procedure get_neighbour_volume
    type(cell_locator) :: nb_cell_location

    call get_neighbour_cell_locator(neighbour_location, nb_cell_location)
    call get_cell_volume(nb_cell_location, V)
  end procedure get_neighbour_volume

  module procedure get_cell_global_index
    associate(mesh => cell_location%mesh)
      if (mesh%nlocal > 0) then ! XXX: Potentially expensive...
        associate(cell => cell_location%cell_idx)
          idxg = mesh%idx_global(cell)
        end associate
      else
        idxg = -1 ! XXX: What should we do in case of too many processors for a given mesh?
      end if
    end associate
  end procedure get_cell_global_index

  module procedure get_neighbour_global_index
    type(cell_locator) :: nb_cell_location
    call get_neighbour_cell_locator(neighbour_location, nb_cell_location)
    call get_global_index(nb_cell_location, nbidxg)
  end procedure get_neighbour_global_index

  module procedure cell_count_neighbours
    associate(mesh => cell_location%mesh, &
         cell => cell_location%cell_idx)
      nnb = mesh%nnb(cell)
    end associate
  end procedure cell_count_neighbours

  module procedure get_neighbour_boundary_status
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
  end procedure get_neighbour_boundary_status

  module procedure get_face_boundary_status
    type(cell_locator) :: cell_location
    type(neighbour_locator) :: neighbour_location

    associate(mesh => face_location%mesh, &
         i => face_location%cell_idx, &
         j => face_location%cell_face_ctr)
      call set_cell_location(mesh, i, cell_location)
      call set_neighbour_location(cell_location, j, neighbour_location)
    end associate
    call get_neighbour_boundary_status(neighbour_location, is_boundary)
  end procedure get_face_boundary_status
  
  
  module procedure get_local_status
    integer :: nbidx

    call get_neighbour_local_index(neighbour_location, nbidx)
    associate(mesh => neighbour_location%mesh)
      if ((nbidx > 0) .and. (nbidx <= mesh%nlocal)) then
        is_local = .true.
      else
        is_local = .false.
      end if
    end associate
  end procedure get_local_status

  module procedure get_cell_local_index
    idx = cell_location%cell_idx
  end procedure get_cell_local_index

  module procedure get_face_local_index
    associate(mesh => face_location%mesh, &
      i => face_location%cell_idx, &
      j => face_location%cell_face_ctr)
      idx = mesh%faceidx(j,i)
    end associate
  end procedure get_face_local_index

  module procedure get_neighbour_local_index
    associate(mesh => neighbour_location%mesh, &
         i => neighbour_location%cell_idx, &
         j => neighbour_location%cell_neighbour_ctr)
      nbidx = mesh%nbidx(j, i)
    end associate
  end procedure get_neighbour_local_index

  subroutine get_neighbour_cell_locator(neighbour_location, cell_location)
    type(neighbour_locator), intent(in) :: neighbour_location
    type(cell_locator), intent(out) :: cell_location

    integer(accs_int) :: nbidx
    
    call get_local_index(neighbour_location, nbidx)
    call set_cell_location(neighbour_location%mesh, nbidx, cell_location)
  end subroutine get_neighbour_cell_locator
end submodule meshing_accessors
