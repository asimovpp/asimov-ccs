module mesh_utils

  use constants, only : ndim
  
  use kinds, only: accs_int, accs_real, accs_err
  use types, only: mesh, face_locator, cell_locator, neighbour_locator
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi
  
  implicit none

  private
  public :: build_square_mesh
  
contains

  function build_square_mesh(nps, l, par_env) result(square_mesh)

    class(parallel_environment) :: par_env
    integer(accs_int), intent(in) :: nps
    real(accs_real), intent(in) :: l

    integer(accs_int) :: istart, iend
    integer(accs_int) :: i, ii, ictr
    integer(accs_int) :: fctr
    integer(accs_int) :: comm_rank, comm_size

    integer(accs_int) :: nbidx, nbidxg

    type(mesh) :: square_mesh

    select type(par_env)
      type is (parallel_environment_mpi)

        square_mesh%nglobal = nps**2            !> (global) Number of cells
        square_mesh%h = l / real(nps, accs_real)

        associate(nglobal=>square_mesh%nglobal, &
                  h=>square_mesh%h)
          
          !! Setup ownership range
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
            iend = iend + 1
          end if

          istart = istart + 1 ! Fortran - 1 indexed
          square_mesh%nlocal = (iend - (istart - 1))

          allocate(square_mesh%idx_global(square_mesh%nlocal))
          allocate(square_mesh%nnb(square_mesh%nlocal))
          allocate(square_mesh%nbidx(4, square_mesh%nlocal))
          allocate(square_mesh%xc(ndim, square_mesh%nlocal))    
          allocate(square_mesh%xf(ndim, 4, square_mesh%nlocal)) !> @note Currently hardcoded as a 2D mesh!
          allocate(square_mesh%vol(square_mesh%nlocal))
          allocate(square_mesh%Af(4, square_mesh%nlocal))    
          allocate(square_mesh%nf(ndim, 4, square_mesh%nlocal)) !> @note Currently hardcoded as a 2D mesh!
          
          square_mesh%nnb(:) = 4 ! All cells have 4 neighbours (possibly ghost/boundary cells)
          square_mesh%vol(:) = square_mesh%h**2 !> @note Mesh is square and 2D
          square_mesh%Af(:, :) = square_mesh%h  !> @note Mesh is square and 2D
          square_mesh%nf(:, :, :) = 0.0_accs_real
          square_mesh%xc(:, :) = 0.0_accs_real
          square_mesh%xf(:, :, :) = 0.0_accs_real
          
          !! Get neighbour indices
          !! XXX: These are global indices and thus may be off-process
          ictr = 1
          do i = istart, iend
            square_mesh%idx_global(ictr) = i
            ii = i - 1

            associate(xc => square_mesh%xc(:, ictr), &
                 xf => square_mesh%xf(:, :, ictr), &
                 nrm => square_mesh%nf(:, :, ictr))
              !! Set cell centre
              xc(1) = (modulo(ii, nps) + 0.5_accs_real) * h
              xc(2) = (ii / nps + 0.5_accs_real) * h
              
              !! Left neighbour
              fctr = 1
              if (modulo(ii, nps) == 0) then
                nbidx = -1
                nbidxg = -1
              else
                nbidx = ictr - 1
                nbidxg = i - 1
              end if
              call build_mesh_add_neighbour(square_mesh, ictr, fctr, nbidx, nbidxg)
              xf(1, fctr) = xc(1) - 0.5_accs_real * h
              xf(2, fctr) = xc(2)
              nrm(1, fctr) = -1.0_accs_real
              nrm(2, fctr) = 0.0_accs_real

              !! Right neighbour
              fctr = 2
              if (modulo(ii, nps) == (nps - 1)) then
                nbidx = -2
                nbidxg = -2
              else
                nbidx = ictr + 1
                nbidxg = i + 1
              end if
              call build_mesh_add_neighbour(square_mesh, ictr, fctr, nbidx, nbidxg)
              xf(1, fctr) = xc(1) + 0.5_accs_real * h
              xf(2, fctr) = xc(2)
              nrm(1, fctr) = 1.0_accs_real
              nrm(2, fctr) = 0.0_accs_real

              !! Down neighbour
              fctr = 3
              if ((ii / nps) == 0) then
                nbidx = -3
                nbidxg = -3
              else
                nbidx = ictr - nps
                nbidxg = i - nps
              end if
              call build_mesh_add_neighbour(square_mesh, ictr, fctr, nbidx, nbidxg)
              xf(1, fctr) = xc(1)
              xf(2, fctr) = xc(2) - 0.5_accs_real * h
              nrm(1, fctr) = 0.0_accs_real
              nrm(2, fctr) = -1.0_accs_real

              !! Up neighbour
              fctr = 4
              if ((ii / nps) == (nps - 1)) then
                nbidx = -4
                nbidxg = -4
              else
                nbidx = ictr + nps
                nbidxg = i + nps
              end if
              call build_mesh_add_neighbour(square_mesh, ictr, fctr, nbidx, nbidxg)
              xf(1, fctr) = xc(1)
              xf(2, fctr) = xc(2) + 0.5_accs_real * h
              nrm(1, fctr) = 0.0_accs_real
              nrm(2, fctr) = 1.0_accs_real
            end associate

            ictr = ictr + 1
          end do
        end associate

      class default
        print *, "Unknown parallel environment type!"
        stop

    end select    
  end function build_square_mesh

  subroutine build_mesh_add_neighbour(meshobj, cellidx, nbctr, nbidx, nbidxg)

    type(mesh), intent(inout) :: meshobj
    integer(accs_int), intent(in) :: cellidx
    integer(accs_int), intent(in) :: nbctr
    integer(accs_int), intent(in) :: nbidx
    integer(accs_int), intent(in) :: nbidxg

    integer(accs_int), dimension(:), allocatable :: tmpidx
    integer(accs_int) :: ng
    logical :: found
    integer(accs_int) :: i
    
    if ((nbidx >= 1) .and. (nbidx <= meshobj%nlocal)) then
      meshobj%nbidx(nbctr, cellidx) = nbidx
    else if (nbidxg < 0) then
      ! Boundary "neighbour" - local index should also be -ve
      if (.not. (nbidx < 0)) then
        print *, "ERROR: boundary neighbours should have -ve indices!"
        stop
      end if
      meshobj%nbidx(nbctr, cellidx) = nbidx
    else
      ng = size(meshobj%idx_global)
      found = .false.
      do i = meshobj%nlocal, ng
        if (meshobj%idx_global(i) == nbidxg) then
          found = .true.
          meshobj%nbidx(nbctr, cellidx) = i
          exit
        end if
      end do

      if (.not. found) then
        allocate(tmpidx(ng + 1))
        tmpidx(1:ng) = meshobj%idx_global(1:ng)
        ng = ng + 1
        tmpidx(ng) = nbidxg
        meshobj%nbidx(nbctr, cellidx) = ng
        
        deallocate(meshobj%idx_global)
        allocate(meshobj%idx_global(ng))
        meshobj%idx_global(:) = tmpidx(:)
        deallocate(tmpidx)
      end if
    end if
    
  end subroutine build_mesh_add_neighbour
  
end module mesh_utils
