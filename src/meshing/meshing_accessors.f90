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
    if (cell_idx > size(geometry%idx_global)) then
      print *, "ERROR: trying to access cell I don't own!", cell_idx, geometry%nlocal
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

  module procedure face_normal
    associate(mesh => face_location%mesh, &
         cell => face_location%cell_idx, &
         face => face_location%cell_face_ctr)
      normal(:) = mesh%nf(:, face, cell)
    end associate
  end procedure face_normal

  module procedure face_area
    associate(mesh => face_location%mesh, &
         cell =>face_location%cell_idx, &
         face =>face_location%cell_face_ctr)
      area = mesh%Af(face, cell)
    end associate
  end procedure face_area

  module procedure cell_centre
    associate(mesh => cell_location%mesh, &
         cell => cell_location%cell_idx)
      x(:) = mesh%xc(:, cell)
    end associate
  end procedure cell_centre

  module procedure face_centre
    associate(mesh => face_location%mesh, &
         cell => face_location%cell_idx, &
         face => face_location%cell_face_ctr)
      x(:) = mesh%xf(:, face, cell)
    end associate
  end procedure face_centre

  module procedure volume
    associate(mesh => cell_location%mesh, &
         cell => cell_location%cell_idx)
      V = mesh%vol(cell)
    end associate
  end procedure volume

  module procedure cell_global_index
    associate(mesh => cell_location%mesh)
      if (mesh%nlocal > 0) then ! XXX: Potentially expensive...
        associate(cell => cell_location%cell_idx)
          idxg = mesh%idx_global(cell)
        end associate
      else
        idxg = -1 ! XXX: What should we do in case of too many processors for a given mesh?
      end if
    end associate
  end procedure cell_global_index

  module procedure neighbour_global_index
    use meshing, only : set_cell_location
    type(cell_locator) :: nb_cell_location
    integer(accs_int) :: nbidx
    call local_index(neighbour_location, nbidx)
    call set_cell_location(nb_cell_location, neighbour_location%mesh, nbidx)
    call global_index(nb_cell_location, nbidxg)
  end procedure neighbour_global_index

  module procedure cell_count_neighbours
    associate(mesh => cell_location%mesh, &
         cell => cell_location%cell_idx)
      nnb = mesh%nnb(cell)
    end associate
  end procedure cell_count_neighbours

  module procedure boundary_status
    integer :: nbidx

    call neighbour_local_index(neighbour_location, nbidx)

    if (nbidx > 0) then
      is_boundary = .false.
    else if (nbidx < 0) then
      is_boundary = .true.
    else
      print *, "ERROR: neighbour index (0) is invalid"
      stop
    end if
  end procedure boundary_status

  module procedure local_status
    integer :: nbidx

    call neighbour_local_index(neighbour_location, nbidx)
    associate(mesh => neighbour_location%mesh)
      if ((nbidx > 0) .and. (nbidx <= mesh%nlocal)) then
        is_local = .true.
      else
        is_local = .false.
      end if
    end associate
  end procedure local_status

  module procedure cell_local_index
    idx = cell_location%cell_idx
  end procedure cell_local_index

  module procedure neighbour_local_index
    associate(mesh => neighbour_location%mesh, &
         i => neighbour_location%cell_idx, &
         j => neighbour_location%cell_neighbour_ctr)
      nbidx = mesh%nbidx(j, i)
    end associate
  end procedure neighbour_local_index

end submodule meshing_accessors
