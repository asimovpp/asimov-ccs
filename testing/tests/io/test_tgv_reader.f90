program test_tgv_reader

  use mpi
  use adios2
  
  implicit none

  ! adios2 handlers
  type(adios2_adios):: adios
  type(adios2_io):: io
  type(adios2_variable):: cell_x
  type(adios2_variable):: fc1
  type(adios2_variable):: fc2  
  type(adios2_engine):: engine
  type(adios2_attribute) :: ncel
  type(adios2_attribute) :: maxfaces
  type(adios2_attribute) :: nfac

  integer(kind=8), dimension(2) :: xyz_sel_start, xyz_sel_count
  integer(kind=8), dimension(1) :: fc_sel_start, fc_sel_count

  real, dimension(:,:), allocatable :: xyz_coords

  integer, dimension(:), allocatable :: face_cell1, face_cell2
  integer, dimension(:), allocatable :: vtxdist   ! Array that indicates rank of vertices local to the process
  integer, dimension(:), allocatable :: xadj      ! Array that points to where local vertices’ adjacency lists begin/end 
  integer, dimension(:), allocatable :: globalbnd  ! Global boundaries
  integer, dimension(:), allocatable :: wkint1d   ! 1D integer work array
  integer, dimension(:,:), allocatable :: wkint2d ! 2D integer work array

  integer :: facenb(2)

  integer :: local_idx_start, local_idx_end, idx_local
  integer :: i, j, k, n, irank, isize, ierr
  integer :: dims
  integer :: num_cells ! total number of grid cells
  integer :: max_faces
  integer :: num_faces
  
  ! Launch MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, isize, ierr)

  ! Number of coordinate dimensions
  dims = 3

  ! Create adios handler passing the communicator, debug mode and error flag
  call adios2_init(adios, MPI_COMM_WORLD, adios2_debug_mode_on, ierr)

  ! Declare an IO process configuration inside adios
  call adios2_declare_io(io, adios, "test_reader", ierr)

  ! Set engine type to HDF5
  call adios2_set_engine(io, "HDF5", ierr)

  ! Open in read mode, this launches an engine
  call adios2_open(engine, io, "TaylorGreen.h5", adios2_mode_read, ierr)

  ! If the file can be opened, proceed with reading contents
  if( ierr == adios2_found ) then

     ! Get the "number of cells" attribute
     call adios2_inquire_attribute(ncel, io, "ncel", ierr)
     call adios2_attribute_data(num_cells, ncel, ierr)
     print*, "Number of cells is: ", num_cells

     ! Store cell range assigned to each process      
     allocate(vtxdist(isize+1))

     vtxdist(1) = 1
     vtxdist(isize + 1) = num_cells + 1
     
     k = int(real(num_cells) / isize)
     j = 1
     do i = 1, isize
        vtxdist(i) = j
        j = j + k
     enddo

     print*, "vtxdist: ", vtxdist

     ! First and last cell index assigned to this process
     local_idx_start = vtxdist(irank + 1)
     local_idx_end = vtxdist(irank + 2) - 1

     print*, "Rank ",irank,", local start and end indices: ", local_idx_start," - ", local_idx_end

     ! Get the xyz coordinates
     call adios2_inquire_variable(cell_x, io, "/cell/x", ierr)

     ! Starting point for reading chunk of data
     xyz_sel_start = (/ 0, vtxdist(irank + 1) - 1 /)
     ! How many data points will be read?
     xyz_sel_count = (/ dims, vtxdist(irank + 2) - vtxdist(irank + 1)/)

     ! Allocate space for chunk of coordinates data of this process
     allocate(xyz_coords(xyz_sel_count(1), xyz_sel_count(2)))
     allocate(xadj(xyz_sel_count(2) + 1))

     ! Select the variable
     call adios2_set_selection(cell_x, 2, xyz_sel_start, xyz_sel_count, ierr )

     ! Read the variable and put data into "xyz_coords"
     call adios2_get(engine, cell_x, xyz_coords, adios2_mode_sync, ierr)

     ! Loop over cells - for debugging only
!     do j=1,xyz_sel_count(2)
!        ! Loop over xyz
!        do i=1,xyz_sel_count(1)
!           if(irank==0) write(*,'(f16.12)',advance='no') xyz_coords(i,j)
!        end do
!        write(*,*) ''
!     end do
     
     ! Get the maximum number of faces
     call adios2_inquire_attribute(maxfaces, io, "maxfaces", ierr)
     call adios2_attribute_data(max_faces, maxfaces, ierr)
     
     print*, "Maximum number of faces is: ", max_faces

     ! Get the number of faces
     call adios2_inquire_attribute(nfac, io, "nfac", ierr)
     call adios2_attribute_data(num_faces, nfac, ierr)
     
     print*, "Number of faces is: ", num_faces

     ! Allocate 2D integer work array
     allocate(wkint2d(vtxdist(irank+2)-vtxdist(irank+1),max_faces+1))
     wkint2d = 0

     ! Allocate global boundaries array
     allocate(globalbnd(num_cells))
     globalbnd = 0

     ! Allocate 1D integer work array
     allocate(wkint1d(isize+1))
     
     wkint1d(1) = 1
     wkint1d(isize+1) = num_faces + 1
     do i=1,isize
        wkint1d(i)=1+(i-1)*int(real(num_faces)/isize)
     enddo
     
!     print*, "wkint1d:", wkint1d

     ! Read /face/cell1 and /face/cell2 variables
     call adios2_inquire_variable(fc1, io, "/face/cell1", ierr)
     call adios2_inquire_variable(fc2, io, "/face/cell2", ierr)

     ! Work out where to start reading data and how many data points
     fc_sel_start = (/ wkint1d(irank+1) - 1 /)
     fc_sel_count = (/ wkint1d(irank+2) - wkint1d(irank+1)/)

!     print*, "Rank ", irank,", fc start: ", fc_sel_start(1)
!     print*, "Rank ", irank,", fc count: ", fc_sel_count(1)

     if(irank == 0) then
        print*,"Building Topology"
     end if

     ! Allocate space for face cell arrays
     allocate(face_cell1(fc_sel_count(1)))
     allocate(face_cell2(fc_sel_count(1)))

     call adios2_set_selection(fc1, 1, fc_sel_start, fc_sel_count, ierr )
     call adios2_set_selection(fc2, 1, fc_sel_start, fc_sel_count, ierr )

     call adios2_get(engine, fc1, face_cell1, adios2_mode_sync, ierr)
     call adios2_get(engine, fc2, face_cell2, adios2_mode_sync, ierr)

!     if(irank==3) print*, face_cell1
!     if(irank==3) print*, face_cell2

     ! Sort through the data
     do n=1,wkint1d(irank+2) - WKint1d(irank+1)
        
        facenb(1) = face_cell1(n)
        facenb(2) = face_cell2(n)
        
        if (facenb(1) .ge. local_idx_start .and. facenb(1) .le. local_idx_end .and. facenb(2) .ne. 0) then
           idx_local = facenb(1) - local_idx_start + 1  ! Local cell index
           k = wkint2d(idx_local, max_faces + 1) + 1    ! Increment number of faces for this cell
           wkint2d(idx_local, k) = facenb(2)            ! Store global index of neighbour cell
           wkint2d(idx_local, max_faces + 1) = k        ! Store number of faces for this cell
        endif

        if (facenb(2) .ge. local_idx_start .and. facenb(2) .le. local_idx_end .and. facenb(1) .ne. 0) then
           idx_local = facenb(2) - local_idx_start + 1  ! Local cell index
           k = wkint2d(idx_local,max_faces + 1) + 1     ! Increment number of faces for this cell
           wkint2d(idx_local, k) = facenb(1)            ! Store global index of neighbour cell
           wkint2d(idx_local, max_faces + 1) = k        ! Store number of faces for this cell
        endif

        ! Boundary face                                                                                                                                                    
        if (facenb(2) .eq. 0) then
           globalbnd(facenb(1)) = globalbnd(facenb(1)) + 1
        end if
        
     enddo

!     print*, "wkint2d: ", wkint2d
!     print*, "globalbnd: ", globalbnd

     if(allocated(xyz_coords)) deallocate(xyz_coords)
     if(allocated(face_cell1)) deallocate(face_cell1)
     if(allocated(face_cell2)) deallocate(face_cell2)
     if(allocated(wkint1d)) deallocate(wkint1d)

     !!!!
     ! NEXT: BUILD PARTITION
     !!!!


     
     if(allocated(vtxdist)) deallocate(vtxdist)
     if(allocated(xadj)) deallocate(xadj)
     if(allocated(wkint2d)) deallocate(wkint2d)
     if(allocated(globalbnd)) deallocate(globalbnd)
  end if

  ! Closes engine1 and deallocates it, becomes unreachable
  call adios2_close(engine, ierr)

  call adios2_finalize(adios, ierr)

  call MPI_Finalize(ierr)

end program test_tgv_reader
