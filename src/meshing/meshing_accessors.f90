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
  
end submodule meshing_accessors
