module mesh_utils

  use constants, only : ndim
  
  use kinds, only: accs_int, accs_real, accs_err
  use types, only: mesh, face_locator, cell_locator, neighbour_locator
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi
  
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
  integer, parameter :: left = 1_accs_int
  integer, parameter :: right = 2_accs_int
  integer, parameter :: down = 3_accs_int
  integer, parameter :: up = 4_accs_int
  
  private
  public :: build_square_mesh
  public :: global_start
  public :: local_count
  
contains

  !> @brief Utility constructor to build a square mesh.
  !
  !> @description Builds a Cartesian grid of NxN cells on the domain LxL.
  !
  !> @param[in] integer(accs_int)    nps         - Number of cells per side of the mesh.
  !> @param[in] real(accs_real)      l           - The length of each side
  !> @param[in] parallel_environment par_env     - The parallel environment to construct the mesh.
  !
  !> @returns   mesh                 square_mesh - The mesh
  function build_square_mesh(nps, l, par_env) result(square_mesh)

    class(parallel_environment) :: par_env
    integer(accs_int), intent(in) :: nps
    real(accs_real), intent(in) :: l

    type(mesh) :: square_mesh

    integer(accs_int) :: istart    !> The (global) starting index of a partition
    integer(accs_int) :: iend      !> The (global) last index of a partition
    integer(accs_int) :: i         !> Loop counter
    integer(accs_int) :: ii        !> Zero-indexed loop counter (simplifies some operations)
    integer(accs_int) :: ictr      !> Local index counter
    integer(accs_int) :: fctr      !> Cell-local face counter
    integer(accs_int) :: comm_rank !> The process ID within the parallel environment
    integer(accs_int) :: comm_size !> The size of the parallel environment

    integer(accs_int) :: nbidx     !> The local index of a neighbour cell
    integer(accs_int) :: nbidxg    !> The global index of a neighbour cell

    select type(par_env)
      type is (parallel_environment_mpi)

        ! Set the global mesh parameters
        square_mesh%nglobal = nps**2            
        square_mesh%h = l / real(nps, accs_real)

        ! Associate aliases to make code easier to read
        associate(nglobal=>square_mesh%nglobal, &
                  h=>square_mesh%h)
          
          ! Determine ownership range
          comm_rank = par_env%proc_id
          comm_size = par_env%num_procs
          istart = global_start(nglobal, par_env%proc_id, par_env%num_procs)
          square_mesh%nlocal = local_count(nglobal, par_env%proc_id, par_env%num_procs)
          iend = istart + (square_mesh%nlocal - 1)

          ! Allocate mesh arrays
          allocate(square_mesh%idx_global(square_mesh%nlocal))
          allocate(square_mesh%nnb(square_mesh%nlocal))
          allocate(square_mesh%nbidx(4, square_mesh%nlocal))
          allocate(square_mesh%xc(ndim, square_mesh%nlocal))    
          allocate(square_mesh%xf(ndim, 4, square_mesh%nlocal)) !> @note Currently hardcoded as a 2D mesh!
          allocate(square_mesh%vol(square_mesh%nlocal))
          allocate(square_mesh%Af(4, square_mesh%nlocal))    
          allocate(square_mesh%nf(ndim, 4, square_mesh%nlocal)) !> @note Currently hardcoded as a 2D mesh!

          ! Initialise mesh arrays
          square_mesh%nnb(:) = 4_accs_int ! All cells have 4 neighbours (possibly ghost/boundary cells)
          square_mesh%vol(:) = square_mesh%h**2 !> @note Mesh is square and 2D
          square_mesh%Af(:, :) = square_mesh%h  !> @note Mesh is square and 2D
          square_mesh%nf(:, :, :) = 0.0_accs_real
          square_mesh%xc(:, :) = 0.0_accs_real
          square_mesh%xf(:, :, :) = 0.0_accs_real

          ! First set the global index of local cells
          ictr = 1_accs_int
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
          ictr = 1_accs_int ! Set local indexing starting from 1...n
          do i = istart, iend 
            ii = i - 1_accs_int

            ! Create aliases for
            ! - xc (centre of cell i)
            ! - xf (centres of faces of cell i)
            ! - nrm (normals of faces of cell i)
            associate(xc => square_mesh%xc(:, ictr), &
                 xf => square_mesh%xf(:, :, ictr), &
                 nrm => square_mesh%nf(:, :, ictr))

              ! Set cell centre
              xc(1) = (modulo(ii, nps) + 0.5_accs_real) * h
              xc(2) = (ii / nps + 0.5_accs_real) * h

              ! Construct left (1) face/neighbour
              fctr = left
              if (modulo(ii, nps) == 0_accs_int) then
                nbidx = -left
                nbidxg = -left
              else
                nbidx = ictr - 1_accs_int
                nbidxg = i - 1_accs_int
              end if
              call build_local_mesh_add_neighbour(square_mesh, ictr, fctr, nbidx, nbidxg)
              xf(1, fctr) = xc(1) - 0.5_accs_real * h
              xf(2, fctr) = xc(2)
              nrm(1, fctr) = -1.0_accs_real
              nrm(2, fctr) = 0.0_accs_real

              ! Construct right (2) face/neighbour
              fctr = right
              if (modulo(ii, nps) == (nps - 1_accs_int)) then
                nbidx = -right
                nbidxg = -right
              else
                nbidx = ictr + 1_accs_int
                nbidxg = i + 1_accs_int
              end if
              call build_local_mesh_add_neighbour(square_mesh, ictr, fctr, nbidx, nbidxg)
              xf(1, fctr) = xc(1) + 0.5_accs_real * h
              xf(2, fctr) = xc(2)
              nrm(1, fctr) = 1.0_accs_real
              nrm(2, fctr) = 0.0_accs_real

              ! Construct down (3) face/neighbour
              fctr = down
              if ((ii / nps) == 0_accs_int) then
                nbidx = -down
                nbidxg = -down
              else
                nbidx = ictr - nps
                nbidxg = i - nps
              end if
              call build_local_mesh_add_neighbour(square_mesh, ictr, fctr, nbidx, nbidxg)
              xf(1, fctr) = xc(1)
              xf(2, fctr) = xc(2) - 0.5_accs_real * h
              nrm(1, fctr) = 0.0_accs_real
              nrm(2, fctr) = -1.0_accs_real

              ! Construct up (4) face/neighbour
              fctr = up
              if ((ii / nps) == (nps - 1_accs_int)) then
                nbidx = -up
                nbidxg = -up
              else
                nbidx = ictr + nps
                nbidxg = i + nps
              end if
              call build_local_mesh_add_neighbour(square_mesh, ictr, fctr, nbidx, nbidxg)
              xf(1, fctr) = xc(1)
              xf(2, fctr) = xc(2) + 0.5_accs_real * h
              nrm(1, fctr) = 0.0_accs_real
              nrm(2, fctr) = 1.0_accs_real
            end associate

            ictr = ictr + 1_accs_int
          end do
        end associate

        square_mesh%ntotal = size(square_mesh%idx_global)
        square_mesh%nhalo = square_mesh%ntotal - square_mesh%nlocal

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
  !> @param[inout] mesh meshobj - the mesh we are assembling neighbours on
  !> @param[in]    integer(accs_int) cellidx - the index of the cell whose neighbours we are assembling
  !> @param[in]    integer(accs_int) nbctr   - the cell-relative neighbour index
  !> @param[in]    integer(accs_int) nbidx   - the local index of the neighbour cell
  !> @param[in]    integer(accs_int) nbidxg  - the global index of the neighbour cell
  subroutine build_local_mesh_add_neighbour(meshobj, cellidx, nbctr, nbidx, nbidxg)

    type(mesh), intent(inout) :: meshobj
    integer(accs_int), intent(in) :: cellidx
    integer(accs_int), intent(in) :: nbctr
    integer(accs_int), intent(in) :: nbidx
    integer(accs_int), intent(in) :: nbidxg

    integer(accs_int) :: ng !> The current number of cells (total = local + halos)
    logical :: found        !> Indicates whether a halo cell was already present
    integer(accs_int) :: i  !> Cell iteration counter
    
    if ((nbidx >= 1_accs_int) .and. (nbidx <= meshobj%nlocal)) then
      ! Neighbour is local
      meshobj%nbidx(nbctr, cellidx) = nbidx
    else if (nbidxg < 0_accs_int) then
      ! Boundary "neighbour" - local index should also be -ve
      if (.not. (nbidx < 0_accs_int)) then
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

    integer(accs_int), intent(in) :: i
    integer(accs_int), dimension(:), allocatable, intent(inout) :: arr ! XXX: Allocatable here be
                                                                       !      dragons! If this were
                                                                       !      intent(out) it would
                                                                       !      be deallocated on entry!

    integer(accs_int) :: n
    integer(accs_int), dimension(:), allocatable :: tmp

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
  
  integer function global_start(n, procid, nproc)

    integer(accs_int), intent(in) :: n
    integer(accs_int), intent(in) :: procid
    integer(accs_int), intent(in) :: nproc

    !! Each PE gets an equal split of the problem with any remainder split equally between the lower
    !! PEs.
    global_start = procid * (n / nproc) + min(procid, modulo(n, nproc))

    !! Fortran indexing
    global_start = global_start + 1
    
  end function global_start

  integer function local_count(n, procid, nproc)

    integer(accs_int), intent(in) :: n
    integer(accs_int), intent(in) :: procid
    integer(accs_int), intent(in) :: nproc

    if (procid < n) then
      local_count = global_start(n, procid, nproc)
      if (procid < (nproc - 1)) then
        local_count = global_start(n, procid + 1, nproc) - local_count
      else
        local_count = n - (local_count - 1)
      end if
    else
      local_count = 0
    end if
    
  end function local_count
  
end module mesh_utils
