submodule (meshing) meshing_accessors

  implicit none
  
contains
  
  module procedure set_face_location
    face_location%mesh => geometry
    face_location%cell_idx = cell_idx
    face_location%cell_face_ctr = cell_face_ctr
  end procedure set_face_location

  module procedure set_cell_location
    ! XXX: Potentially expensive...
    if (cell_idx > geometry%ntotal) then
      print *, "ERROR: trying to access cell I don't have access to!", cell_idx, geometry%nlocal
      stop
    else
      cell_location%mesh => geometry
      cell_location%cell_idx = cell_idx
    end if
  end procedure set_cell_location
  
  module procedure set_neighbour_location
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
  end procedure set_neighbour_location

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

  module procedure get_face_centre
    associate(mesh => face_location%mesh, &
         cell => face_location%cell_idx, &
         face => face_location%cell_face_ctr)
      x(:) = mesh%xf(:, face, cell)
    end associate
  end procedure get_face_centre

  module procedure get_volume
    associate(mesh => cell_location%mesh, &
         cell => cell_location%cell_idx)
      V = mesh%vol(cell)
    end associate
  end procedure get_volume

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
    integer(accs_int) :: nbidx
    call get_local_index(neighbour_location, nbidx)
    call set_cell_location(nb_cell_location, neighbour_location%mesh, nbidx)
    call get_global_index(nb_cell_location, nbidxg)
  end procedure get_neighbour_global_index

  module procedure cell_count_neighbours
    associate(mesh => cell_location%mesh, &
         cell => cell_location%cell_idx)
      nnb = mesh%nnb(cell)
    end associate
  end procedure cell_count_neighbours

  module procedure get_boundary_status
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
  end procedure get_boundary_status

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

  module procedure get_neighbour_local_index
    associate(mesh => neighbour_location%mesh, &
         i => neighbour_location%cell_idx, &
         j => neighbour_location%cell_neighbour_ctr)
      nbidx = mesh%nbidx(j, i)
    end associate
  end procedure get_neighbour_local_index

end submodule meshing_accessors
