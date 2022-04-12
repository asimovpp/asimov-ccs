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

  !> @note Named constants for faces of hexahedral cells follow the convention that the lower
  !!       boundary on a given axis is numbered first, i.e.
  !!
  !!           4
  !!     +----------+
  !!     |          |
  !!     |          |
  !!   1 |          | 2
  !!     |          |
  !!     +----------+
  !!           3
  !!
  integer, parameter :: left = 1_ccs_int
  integer, parameter :: right = 2_ccs_int
  integer, parameter :: down = 3_ccs_int
  integer, parameter :: up = 4_ccs_int
  
  private
  public :: build_square_mesh
  public :: count_mesh_faces
  
contains

  !> @brief Utility constructor to build a square mesh.
  !
  !> @description Builds a Cartesian grid of NxN cells on the domain LxL.
  !
  !> @param[in] integer(ccs_int)    nps         - Number of cells per side of the mesh.
  !> @param[in] real(ccs_real)      l           - The length of each side
  !> @param[in] parallel_environment par_env     - The parallel environment to construct the mesh.
  !
  !> @returns   mesh                 square_mesh - The mesh
  function build_square_mesh(par_env, nps, l) result(square_mesh)

    class(parallel_environment), intent(in) :: par_env
    integer(ccs_int), intent(in) :: nps
    real(ccs_real), intent(in) :: l

    type(ccs_mesh) :: square_mesh

    integer(ccs_int) :: istart    !> The (global) starting index of a partition
    integer(ccs_int) :: iend      !> The (global) last index of a partition
    integer(ccs_int) :: i         !> Loop counter
    integer(ccs_int) :: ii        !> Zero-indexed loop counter (simplifies some operations)
    integer(ccs_int) :: ictr      !> Local index counter
    integer(ccs_int) :: fctr      !> Cell-local face counter
    integer(ccs_int) :: comm_rank !> The process ID within the parallel environment
    integer(ccs_int) :: comm_size !> The size of the parallel environment

    integer(ccs_int) :: nbidx     !> The local index of a neighbour cell
    integer(ccs_int) :: nbidxg    !> The global index of a neighbour cell

    select type(par_env)
      type is (parallel_environment_mpi)

        ! Set the global mesh parameters
        square_mesh%nglobal = nps**2            
        square_mesh%h = l / real(nps, ccs_real)

        ! Associate aliases to make code easier to read
        associate(nglobal=>square_mesh%nglobal, &
                  h=>square_mesh%h)
          
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
          square_mesh%nlocal = (iend - (istart - 1_ccs_int))

          ! Allocate mesh arrays
          allocate(square_mesh%idx_global(square_mesh%nlocal))
          allocate(square_mesh%nnb(square_mesh%nlocal))
          allocate(square_mesh%nbidx(4, square_mesh%nlocal))
          allocate(square_mesh%faceidx(4, square_mesh%nlocal))

          ! Initialise mesh arrays
          square_mesh%nnb(:) = 4_ccs_int ! All cells have 4 neighbours (possibly ghost/boundary cells)

          ! First set the global index of local cells
          ictr = 1_ccs_int
          do i = istart, iend
            square_mesh%idx_global(ictr) = i
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
              nbidx = -left
              nbidxg = -left
            else
              nbidx = ictr - 1_ccs_int
              nbidxg = i - 1_ccs_int
            end if
            call build_local_mesh_add_neighbour(ictr, fctr, nbidx, nbidxg, square_mesh)

            ! Construct right (2) face/neighbour
            fctr = right
            if (modulo(ii, nps) == (nps - 1_ccs_int)) then
              nbidx = -right
              nbidxg = -right
            else
              nbidx = ictr + 1_ccs_int
              nbidxg = i + 1_ccs_int
            end if
            call build_local_mesh_add_neighbour(ictr, fctr, nbidx, nbidxg, square_mesh)

            ! Construct down (3) face/neighbour
            fctr = down
            if ((ii / nps) == 0_ccs_int) then
              nbidx = -down
              nbidxg = -down
            else
              nbidx = ictr - nps
              nbidxg = i - nps
            end if
            call build_local_mesh_add_neighbour(ictr, fctr, nbidx, nbidxg, square_mesh)

            ! Construct up (4) face/neighbour
            fctr = up
            if ((ii / nps) == (nps - 1_ccs_int)) then
              nbidx = -up
              nbidxg = -up
            else
              nbidx = ictr + nps
              nbidxg = i + nps
            end if
            call build_local_mesh_add_neighbour(ictr, fctr, nbidx, nbidxg, square_mesh)

            ictr = ictr + 1_ccs_int
          end do
        end associate

        square_mesh%ntotal = size(square_mesh%idx_global)
        square_mesh%nhalo = square_mesh%ntotal - square_mesh%nlocal

        allocate(square_mesh%xc(ndim, square_mesh%ntotal))    
        allocate(square_mesh%xf(ndim, 4, square_mesh%nlocal)) !> @note Currently hardcoded as a 2D mesh!
        allocate(square_mesh%vol(square_mesh%ntotal))
        allocate(square_mesh%Af(4, square_mesh%nlocal))    
        allocate(square_mesh%nf(ndim, 4, square_mesh%nlocal)) !> @note Currently hardcoded as a 2D mesh!

        square_mesh%vol(:) = square_mesh%h**2 !> @note Mesh is square and 2D
        square_mesh%nf(:, :, :) = 0.0_ccs_real
        square_mesh%xc(:, :) = 0.0_ccs_real
        square_mesh%xf(:, :, :) = 0.0_ccs_real
        square_mesh%Af(:, :) = square_mesh%h  !> @note Mesh is square and 2D

        associate(h => square_mesh%h)
          do i = 1_ccs_int, square_mesh%ntotal
            ii = square_mesh%idx_global(i)

            associate(xc => square_mesh%xc(:, i))
              ! Set cell centre
              xc(1) = (modulo(ii-1, nps) + 0.5_ccs_real) * h
              xc(2) = ((ii - 1) / nps + 0.5_ccs_real) * h
            end associate
          end do

          do i = 1_ccs_int, square_mesh%nlocal
            associate(xc => square_mesh%xc(:, i), &
                 xf => square_mesh%xf(:, :, i), &
                 nrm => square_mesh%nf(:, :, i))

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

        square_mesh%nfaces_local = count_mesh_faces(square_mesh)

        call set_cell_face_indices(square_mesh)

      class default
        print *, "Unknown parallel environment type!"
        stop

    end select    
  end function build_square_mesh

  !> @brief Helper subroutine to add a neighbour to a cell's neighbour list.
  !
  !> @description Given a local and global index for a neighbour there are 3 possibilities:
  !!              1) the local and the neighbour is added immediately
  !!              2) the global index is negative indicating it is a boundary and the "neighbour" is
  !!                 added immediately
  !!              3) the index is not local:
  !!                 a) the global index is already in the off-process list (halos), the neighbour
  !!                    is added immediately
  !!                 b) this is a new halo cell, the list of global indices must be grown to
  !!                    accomodate before adding the neighbour.
  !
  !> @param[in]    integer(ccs_int) cellidx - the index of the cell whose neighbours we are assembling
  !> @param[in]    integer(ccs_int) nbctr   - the cell-relative neighbour index
  !> @param[in]    integer(ccs_int) nbidx   - the local index of the neighbour cell
  !> @param[in]    integer(ccs_int) nbidxg  - the global index of the neighbour cell
  !> @param[inout] mesh meshobj - the mesh we are assembling neighbours on
  subroutine build_local_mesh_add_neighbour(cellidx, nbctr, nbidx, nbidxg, meshobj)

    integer(ccs_int), intent(in) :: cellidx
    integer(ccs_int), intent(in) :: nbctr
    integer(ccs_int), intent(in) :: nbidx
    integer(ccs_int), intent(in) :: nbidxg
    type(ccs_mesh), intent(inout) :: meshobj

    integer(ccs_int) :: ng !> The current number of cells (total = local + halos)
    logical :: found        !> Indicates whether a halo cell was already present
    integer(ccs_int) :: i  !> Cell iteration counter
    
    if ((nbidx >= 1_ccs_int) .and. (nbidx <= meshobj%nlocal)) then
      ! Neighbour is local
      meshobj%nbidx(nbctr, cellidx) = nbidx
    else if (nbidxg < 0_ccs_int) then
      ! Boundary "neighbour" - local index should also be -ve
      if (.not. (nbidx < 0_ccs_int)) then
        print *, "ERROR: boundary neighbours should have -ve indices!"
        stop
      end if
      meshobj%nbidx(nbctr, cellidx) = nbidx
    else
      ! Neighbour is in a halo

      ! First check if neighbour is already present in halo
      ng = size(meshobj%idx_global)
      found = .false.
      do i = meshobj%nlocal + 1, ng
        if (meshobj%idx_global(i) == nbidxg) then
          found = .true.
          meshobj%nbidx(nbctr, cellidx) = i
          exit
        end if
      end do

      ! If neighbour was not present append to global index list (the end of the global index list
      ! becoming its local index).
      ! XXX: Note this currently copies into an n+1 temporary, reallocates and then copies back to
      !      the (extended) original array.
      if (.not. found) then
        if ((ng + 1) > meshobj%nglobal) then
          print *, "ERROR: Trying to create halo that exceeds global mesh size!"
          stop
        end if
        
        call append_to_arr(nbidxg, meshobj%idx_global)
        ng = size(meshobj%idx_global)
        meshobj%nbidx(nbctr, cellidx) = ng
      end if
    end if
    
  end subroutine build_local_mesh_add_neighbour

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

  !> @brief Count the number of faces in the mesh
  !
  !> @param[in]  cell_mesh - the mesh
  !> @param[out] nfaces    - number of cell faces
  function count_mesh_faces(cell_mesh) result(nfaces)

    ! Arguments
    type(ccs_mesh), intent(in) :: cell_mesh

    ! Result
    integer(ccs_int) :: nfaces

    ! Local variables
    type(cell_locator) :: self_loc
    type(neighbour_locator) :: loc_nb
    integer(ccs_int) :: self_idx, local_idx
    integer(ccs_int) :: j
    integer(ccs_int) :: nnb
    integer(ccs_int) :: nfaces_int       !> Internal face count
    integer(ccs_int) :: nfaces_bnd       !> Boundary face count
    integer(ccs_int) :: nfaces_interface !> Process interface face count
    logical :: is_boundary
    logical :: is_local

    ! Initialise
    nfaces_int = 0
    nfaces_bnd = 0
    nfaces_interface = 0

    ! Loop over cells
    do local_idx = 1, cell_mesh%nlocal
      call set_cell_location(cell_mesh, local_idx, self_loc)
      call get_global_index(self_loc, self_idx)
      call count_neighbours(self_loc, nnb)

      do j = 1, nnb
        call set_neighbour_location(self_loc, j, loc_nb)
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

  subroutine set_cell_face_indices(cell_mesh)

    ! Arguments
    type(ccs_mesh), intent(inout) :: cell_mesh

    ! Local variables
    type(cell_locator) :: self_loc !> Current cell
    type(neighbour_locator) :: loc_nb !> Neighbour
    integer(ccs_int) :: index_nb, local_idx
    integer(ccs_int) :: face_idx
    integer(ccs_int) :: nnb
    integer(ccs_int) :: j
    integer(ccs_int) :: icnt  !> Face index counter
    logical :: is_boundary

    icnt = 0

    ! Loop over cells
    do local_idx = 1, cell_mesh%nlocal
      call set_cell_location(cell_mesh, local_idx, self_loc)
      call count_neighbours(self_loc, nnb)

      do j = 1, nnb
        call set_neighbour_location(self_loc, j, loc_nb)
        call get_local_index(loc_nb, index_nb)
        call get_boundary_status(loc_nb, is_boundary)

        if (.not. is_boundary) then
          ! Cell with lowest local index assigns an index to the face
          if (local_idx < index_nb) then
            icnt = icnt + 1
            call set_face_index(local_idx, j, icnt, cell_mesh)
          else
            ! Find corresponding face in neighbour cell
            ! (To be improved, this seems inefficient!)
            face_idx = get_neighbour_face_index(cell_mesh, local_idx, index_nb)
            call set_face_index(local_idx, j, face_idx, cell_mesh)
          endif
        else
          icnt = icnt + 1
          call set_face_index(local_idx, j, icnt, cell_mesh)
        endif
      end do  ! End loop over current cell's neighbours
    end do    ! End loop over local cells

  end subroutine set_cell_face_indices

  !> @brief Computes the index of the face shared by the cells denoted by the specified 
  !!        local index and neighbouring index
  !<
  !> @param[in] mesh      - the mesh
  !> @param[in] local_idx - the current cell index
  !> @param[in] index_nb  - the index of the neighbouring cell
  function get_neighbour_face_index(mesh, local_idx, index_nb) result(face_idx)
    type(ccs_mesh), intent(in) :: mesh
    integer(ccs_int), intent(in) :: local_idx
    integer(ccs_int), intent(in) :: index_nb
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
      if (index_nb_nb == local_idx) then
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
