module mesh_utils

  use constants, only : ndim
  
  use kinds, only: ccs_int, ccs_real
  use types, only: ccs_mesh, face_locator, cell_locator, neighbour_locator
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi
  use meshing, only: get_global_index, get_local_index, count_neighbours, &
                     set_cell_location, set_neighbour_location, set_face_location, &
                     set_face_index, get_boundary_status, get_local_status

  
  implicit none

  !^ @note Named constants for faces of hexahedral cells follow the convention that the lower
  !        boundary on a given axis is numbered first, i.e.
  ! 
  !             4
  !       +----------+
  !       |          |
  !       |          |
  !     1 |          | 2
  !       |          |
  !       +----------+
  !             3
  !
  !  @endnote 
  integer, parameter :: left = 1_ccs_int
  integer, parameter :: right = 2_ccs_int
  integer, parameter :: down = 3_ccs_int
  integer, parameter :: up = 4_ccs_int
  
  private
  public :: build_square_mesh
  public :: count_mesh_faces
  
contains

  !v Utility constructor to build a square mesh.
  !
  !  Builds a Cartesian grid of NxN cells on the domain LxL.
  function build_square_mesh(par_env, nps, l) result(mesh)

    class(parallel_environment), intent(in) :: par_env !< Number of cells per side of the mesh.
    integer(ccs_int), intent(in) :: nps                !< The length of each side.
    real(ccs_real), intent(in) :: l                    !< The parallel environment to construct the mesh.

    type(ccs_mesh) :: mesh !< The resulting mesh.

    integer(ccs_int) :: istart          ! The (global) starting index of a partition
    integer(ccs_int) :: iend            ! The (global) last index of a partition
    integer(ccs_int) :: i               ! Loop counter
    integer(ccs_int) :: ii              ! Zero-indexed loop counter (simplifies some operations)
    integer(ccs_int) :: ictr            ! Local index counter
    integer(ccs_int) :: fctr            ! Cell-local face counter
    integer(ccs_int) :: comm_rank       ! The process ID within the parallel environment
    integer(ccs_int) :: comm_size       ! The size of the parallel environment

    integer(ccs_int) :: index_nb        ! The local index of a neighbour cell
    integer(ccs_int) :: global_index_nb ! The global index of a neighbour cell

    select type(par_env)
      type is (parallel_environment_mpi)

        ! Set the global mesh parameters
        mesh%nglobal = nps**2            
        mesh%h = l / real(nps, ccs_real)

        ! Associate aliases to make code easier to read
        associate(nglobal=>mesh%nglobal, &
                  h=>mesh%h)
          
          ! Determine ownership range (based on PETSc ex3.c)
          comm_rank = par_env%proc_id
          comm_size = par_env%num_procs
          istart = comm_rank * (nglobal / comm_size)
          if (modulo(nglobal, comm_size) < comm_rank) then
            istart = istart + modulo(nglobal, comm_size)
          else
            istart = istart + comm_rank
          end if
          iend = istart + nglobal / comm_size
          if (modulo(nglobal, comm_size) > comm_rank) then
            iend = iend + 1_ccs_int
          end if

          ! Fix indexing and determine size of local partition
          istart = istart + 1_ccs_int
          mesh%nlocal = (iend - (istart - 1_ccs_int))

          ! Allocate mesh arrays
          allocate(mesh%idx_global(mesh%nlocal))
          allocate(mesh%nnb(mesh%nlocal))
          allocate(mesh%index_nb(4, mesh%nlocal))
          allocate(mesh%faceidx(4, mesh%nlocal))

          ! Initialise mesh arrays
          mesh%nnb(:) = 4_ccs_int ! All cells have 4 neighbours (possibly ghost/boundary cells)

          ! First set the global index of local cells
          ictr = 1_ccs_int
          do i = istart, iend
            mesh%idx_global(ictr) = i
            ictr = ictr + 1
          end do
          
          ! Assemble cells and faces
          ! XXX: Negative neighbour indices are used to indicate boundaries using the same numbering
          !      as cell-relative neighbour indexing, i.e.
          !        -1 = left boundary
          !        -2 = right boundary
          !        -3 = down boundary
          !        -4 = up boundary
          ictr = 1_ccs_int ! Set local indexing starting from 1...n
          do i = istart, iend 
            ii = i - 1_ccs_int

            ! Construct left (1) face/neighbour
            fctr = left
            if (modulo(ii, nps) == 0_ccs_int) then
              index_nb = -left
              global_index_nb = -left
            else
              index_nb = ictr - 1_ccs_int
              global_index_nb = i - 1_ccs_int
            end if
            call build_local_mesh_add_neighbour(ictr, fctr, index_nb, global_index_nb, mesh)

            ! Construct right (2) face/neighbour
            fctr = right
            if (modulo(ii, nps) == (nps - 1_ccs_int)) then
              index_nb = -right
              global_index_nb = -right
            else
              index_nb = ictr + 1_ccs_int
              global_index_nb = i + 1_ccs_int
            end if
            call build_local_mesh_add_neighbour(ictr, fctr, index_nb, global_index_nb, mesh)

            ! Construct down (3) face/neighbour
            fctr = down
            if ((ii / nps) == 0_ccs_int) then
              index_nb = -down
              global_index_nb = -down
            else
              index_nb = ictr - nps
              global_index_nb = i - nps
            end if
            call build_local_mesh_add_neighbour(ictr, fctr, index_nb, global_index_nb, mesh)

            ! Construct up (4) face/neighbour
            fctr = up
            if ((ii / nps) == (nps - 1_ccs_int)) then
              index_nb = -up
              global_index_nb = -up
            else
              index_nb = ictr + nps
              global_index_nb = i + nps
            end if
            call build_local_mesh_add_neighbour(ictr, fctr, index_nb, global_index_nb, mesh)

            ictr = ictr + 1_ccs_int
          end do
        end associate

        mesh%ntotal = size(mesh%idx_global)
        mesh%nhalo = mesh%ntotal - mesh%nlocal

        allocate(mesh%xc(ndim, mesh%ntotal))    
        allocate(mesh%xf(ndim, 4, mesh%nlocal)) !< @note Currently hardcoded as a 2D mesh! @endnote
        allocate(mesh%vol(mesh%ntotal))
        allocate(mesh%Af(4, mesh%nlocal))    
        allocate(mesh%nf(ndim, 4, mesh%nlocal)) !        Currently hardcoded as a 2D mesh!

        mesh%vol(:) = mesh%h**2 !< @note Mesh is square and 2D @endnote
        mesh%nf(:, :, :) = 0.0_ccs_real
        mesh%xc(:, :) = 0.0_ccs_real
        mesh%xf(:, :, :) = 0.0_ccs_real
        mesh%Af(:, :) = mesh%h  !        Mesh is square and 2D

        associate(h => mesh%h)
          do i = 1_ccs_int, mesh%ntotal
            ii = mesh%idx_global(i)

            associate(xc => mesh%xc(:, i))
              ! Set cell centre
              xc(1) = (modulo(ii-1, nps) + 0.5_ccs_real) * h
              xc(2) = ((ii - 1) / nps + 0.5_ccs_real) * h
            end associate
          end do

          do i = 1_ccs_int, mesh%nlocal
            associate(xc => mesh%xc(:, i), &
                 xf => mesh%xf(:, :, i), &
                 nrm => mesh%nf(:, :, i))

              fctr = left
              xf(1, fctr) = xc(1) - 0.5_ccs_real * h
              xf(2, fctr) = xc(2)
              nrm(1, fctr) = -1.0_ccs_real
              nrm(2, fctr) = 0.0_ccs_real

              fctr = right
              xf(1, fctr) = xc(1) + 0.5_ccs_real * h
              xf(2, fctr) = xc(2)
              nrm(1, fctr) = 1.0_ccs_real
              nrm(2, fctr) = 0.0_ccs_real
              
              fctr = down
              xf(1, fctr) = xc(1)
              xf(2, fctr) = xc(2) - 0.5_ccs_real * h
              nrm(1, fctr) = 0.0_ccs_real
              nrm(2, fctr) = -1.0_ccs_real

              fctr = up
              xf(1, fctr) = xc(1)
              xf(2, fctr) = xc(2) + 0.5_ccs_real * h
              nrm(1, fctr) = 0.0_ccs_real
              nrm(2, fctr) = 1.0_ccs_real
            end associate
          end do
        end associate

        mesh%nfaces_local = count_mesh_faces(mesh)

        call set_cell_face_indices(mesh)

      class default
        print *, "Unknown parallel environment type!"
        stop

    end select    
  end function build_square_mesh

  !v Helper subroutine to add a neighbour to a cell's neighbour list.
  !
  !  Given a local and global index for a neighbour there are 3 possibilities:
  !  1) the local and the neighbour is added immediately
  !  2) the global index is negative indicating it is a boundary and the "neighbour" is
  !     added immediately
  !  3) the index is not local:
  !     a) the global index is already in the off-process list (halos), the neighbour
  !        is added immediately
  !     b) this is a new halo cell, the list of global indices must be grown to
  !        accomodate before adding the neighbour.
  subroutine build_local_mesh_add_neighbour(index_p, index_p_nb, index_nb, global_index_nb, mesh)

    integer(ccs_int), intent(in) :: index_p !< the index of the cell whose neighbours we are assembling 
    integer(ccs_int), intent(in) :: index_p_nb !< the cell-relative neighbour index
    integer(ccs_int), intent(in) :: index_nb !< the local index of the neighbour cell
    integer(ccs_int), intent(in) :: global_index_nb !< the global index of the neighbour cell
    type(ccs_mesh), intent(inout) :: mesh !< the mesh we are assembling neighbours on

    integer(ccs_int) :: ng  ! The current number of cells (total = local + halos)
    logical :: found        ! Indicates whether a halo cell was already present
    integer(ccs_int) :: i   ! Cell iteration counter
    
    if ((index_nb >= 1_ccs_int) .and. (index_nb <= mesh%nlocal)) then
      ! Neighbour is local
      mesh%index_nb(index_p_nb, index_p) = index_nb
    else if (global_index_nb < 0_ccs_int) then
      ! Boundary "neighbour" - local index should also be -ve
      if (.not. (index_nb < 0_ccs_int)) then
        print *, "ERROR: boundary neighbours should have -ve indices!"
        stop
      end if
      mesh%index_nb(index_p_nb, index_p) = index_nb
    else
      ! Neighbour is in a halo

      ! First check if neighbour is already present in halo
      ng = size(mesh%idx_global)
      found = .false.
      do i = mesh%nlocal + 1, ng
        if (mesh%idx_global(i) == global_index_nb) then
          found = .true.
          mesh%index_nb(index_p_nb, index_p) = i
          exit
        end if
      end do

      ! If neighbour was not present append to global index list (the end of the global index list
      ! becoming its local index).
      ! XXX: Note this currently copies into an n+1 temporary, reallocates and then copies back to
      !      the (extended) original array.
      if (.not. found) then
        if ((ng + 1) > mesh%nglobal) then
          print *, "ERROR: Trying to create halo that exceeds global mesh size!"
          stop
        end if
        
        call append_to_arr(global_index_nb, mesh%idx_global)
        ng = size(mesh%idx_global)
        mesh%index_nb(index_p_nb, index_p) = ng
      end if
    end if
    
  end subroutine build_local_mesh_add_neighbour

  !v @note Docs needed.
  subroutine append_to_arr(i, arr)

    integer(ccs_int), intent(in) :: i
    integer(ccs_int), dimension(:), allocatable, intent(inout) :: arr ! XXX: Allocatable here be
                                                                      !      dragons! If this were
                                                                      !      intent(out) it would
                                                                      !      be deallocated on entry!
    integer(ccs_int) :: n
    integer(ccs_int), dimension(:), allocatable :: tmp

    n = size(arr)

    allocate(tmp(n + 1))

    tmp(1:n) = arr(1:n)

    n = n + 1
    tmp(n) = i

    deallocate(arr)
    allocate(arr(n))
    arr(:) = tmp(:)
    deallocate(tmp)
    
  end subroutine append_to_arr

  !v Count the number of faces in the mesh
  function count_mesh_faces(mesh) result(nfaces)

    ! Arguments
    type(ccs_mesh), intent(in) :: mesh !< the mesh

    ! Result
    integer(ccs_int) :: nfaces !< number of cell faces

    ! Local variables
    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb
    integer(ccs_int) :: global_index_p, index_p
    integer(ccs_int) :: j
    integer(ccs_int) :: nnb
    integer(ccs_int) :: nfaces_int       ! Internal face count
    integer(ccs_int) :: nfaces_bnd       ! Boundary face count
    integer(ccs_int) :: nfaces_interface ! Process interface face count
    logical :: is_boundary
    logical :: is_local

    ! Initialise
    nfaces_int = 0
    nfaces_bnd = 0
    nfaces_interface = 0

    ! Loop over cells
    do index_p = 1, mesh%nlocal
      call set_cell_location(mesh, index_p, loc_p)
      call get_global_index(loc_p, global_index_p)
      call count_neighbours(loc_p, nnb)

      do j = 1, nnb
        call set_neighbour_location(loc_p, j, loc_nb)
        call get_boundary_status(loc_nb, is_boundary)

        if (.not. is_boundary) then
          call get_local_status(loc_nb, is_local)
          
          if (is_local) then
            ! Interior face
            nfaces_int = nfaces_int + 1
          else
            ! Process boundary face
            nfaces_interface = nfaces_interface + 1
          end if
        else
          ! Boundary face
          nfaces_bnd = nfaces_bnd + 1
        endif
      end do
    end do

    ! Interior faces will be counted twice
    nfaces = (nfaces_int / 2) + nfaces_interface + nfaces_bnd

  end function count_mesh_faces

  !v @note Docs needed.
  subroutine set_cell_face_indices(mesh)

    ! Arguments
    type(ccs_mesh), intent(inout) :: mesh

    ! Local variables
    type(cell_locator) :: loc_p           ! Current cell
    type(neighbour_locator) :: loc_nb     ! Neighbour
    integer(ccs_int) :: index_nb, index_p
    integer(ccs_int) :: face_idx
    integer(ccs_int) :: nnb
    integer(ccs_int) :: j
    integer(ccs_int) :: icnt              ! Face index counter
    logical :: is_boundary

    icnt = 0

    ! Loop over cells
    do index_p = 1, mesh%nlocal
      call set_cell_location(mesh, index_p, loc_p)
      call count_neighbours(loc_p, nnb)

      do j = 1, nnb
        call set_neighbour_location(loc_p, j, loc_nb)
        call get_local_index(loc_nb, index_nb)
        call get_boundary_status(loc_nb, is_boundary)

        if (.not. is_boundary) then
          ! Cell with lowest local index assigns an index to the face
          if (index_p < index_nb) then
            icnt = icnt + 1
            call set_face_index(index_p, j, icnt, mesh)
          else
            ! Find corresponding face in neighbour cell
            ! (To be improved, this seems inefficient!)
            face_idx = get_neighbour_face_index(mesh, index_p, index_nb)
            call set_face_index(index_p, j, face_idx, mesh)
          endif
        else
          icnt = icnt + 1
          call set_face_index(index_p, j, icnt, mesh)
        endif
      end do  ! End loop over current cell's neighbours
    end do    ! End loop over local cells

  end subroutine set_cell_face_indices

  !v Computes the index of the face shared by the cells denoted by the specified 
  !  local index and neighbouring index
  function get_neighbour_face_index(mesh, index_p, index_nb) result(face_idx)
    type(ccs_mesh), intent(in) :: mesh       !< the mesh
    integer(ccs_int), intent(in) :: index_p  !< the current cell index
    integer(ccs_int), intent(in) :: index_nb !< the index of the neighbouring cell
    integer(ccs_int) :: face_idx

    ! Local variables
    integer(ccs_int) :: k
    integer(ccs_int) :: nnb_nb
    type(cell_locator) :: loc_nb 
    type(neighbour_locator) :: loc_nb_nb
    type(face_locator) :: face_loc
    integer(ccs_int) :: index_nb_nb

    call set_cell_location(mesh, index_nb, loc_nb)
    call count_neighbours(loc_nb, nnb_nb)
    do k = 1, nnb_nb
      call set_neighbour_location(loc_nb, k, loc_nb_nb)
      call get_local_index(loc_nb_nb, index_nb_nb)
      if (index_nb_nb == index_p) then
        call set_face_location(mesh, index_nb, k, face_loc)
        call get_local_index(face_loc, face_idx)
        exit ! Exit the loop, as found shared face
      else if (k == nnb_nb) then
        print *, "ERROR: Failed to find face in owning cell"
        stop 1
      endif
    end do
  end function get_neighbour_face_index
  
end module mesh_utils
