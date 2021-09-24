module accs_mesh

  use accs_kinds
  use accs_types
  
  implicit none

  private
  public :: build_mesh
  
contains

  function build_mesh(nps, l, comm) result(square_mesh)

    use mpi
    
    integer(accs_int), intent(in) :: nps
    real(accs_real), intent(in) :: l
    integer, intent(in) :: comm

    integer(accs_int) :: istart, iend
    integer(accs_int) :: i, ii, ictr
    type(mesh) :: square_mesh

    integer :: comm_rank, comm_size
    integer(accs_err) :: ierr

    square_mesh%n = nps**2            ! Number of cells
    square_mesh%h = l / nps
    associate(n=>square_mesh%n, &
         h=>square_mesh%h)
      square_mesh%Af = square_mesh%h     ! Face area
      square_mesh%vol = square_mesh%h**2 ! Cell volume

      !! Setup ownership range
      call MPI_Comm_size(comm, comm_size, ierr)
      call MPI_Comm_rank(comm, comm_rank, ierr)
      istart = comm_rank * (n / comm_size)
      if (modulo(square_mesh%n, comm_size) < comm_rank) then
         istart = istart + modulo(n, comm_size)
      else
         istart = istart + comm_rank
      end if
      istart = istart + 1 ! Fortran - 1 indexed
  
      iend = istart + n / comm_size
      if (modulo(square_mesh%n, comm_size) > comm_rank) then
         iend = iend + 1
      end if
      iend = iend - 1

      square_mesh%nlocal = (iend - istart) + 1
      allocate(square_mesh%idx_global(square_mesh%nlocal), square_mesh%nnb(square_mesh%nlocal))
      square_mesh%nnb(:) = 4                    ! All cells have 4 neighbours (possibly ghost/boundary cells)
      allocate(square_mesh%nbidx(4, square_mesh%nlocal))
    
      !! Get neighbour indices
      !! XXX: These are global indices and thus may be off-process
      ictr = 1
      do i = istart, iend
         square_mesh%idx_global(ictr) = i
         ii = i - 1
       
         !! Left neighbour
         if (modulo(ii, nps) == 0) then
            square_mesh%nbidx(1, ictr) = -1
         else
            square_mesh%nbidx(1, ictr) = i - 1
         end if

         !! Right neighbour
         if (modulo(ii, nps) == (nps - 1)) then
            square_mesh%nbidx(2, ictr) = -2
         else
            square_mesh%nbidx(2, ictr) = i + 1
         end if

         !! Down neighbour
         if ((ii / nps) == 0) then
            square_mesh%nbidx(3, ictr) = -3
         else
            square_mesh%nbidx(3, ictr) = i - nps
         end if

         !! Up neighbour
         if ((ii / nps) == (nps - 1)) then
            square_mesh%nbidx(4, ictr) = -4
         else
            square_mesh%nbidx(4, ictr) = i + nps
         end if

         ictr = ictr + 1
      end do
    end associate
    
  end function build_mesh
end module accs_mesh
