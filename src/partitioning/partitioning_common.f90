submodule (partitioning) partitioning_common
#include "ccs_macros.inc"

  use kinds, only: ccs_int
  use utils, only : str, debug_print
  use parallel_types_mpi, only: parallel_environment_mpi

implicit none

  contains

  !v Read the topology data from an input (HDF5) file
  ! This subroutine assumes the following names are used in the file:
  ! "ncel" - the total number of cells
  ! "nfac" - the total number of faces
  ! "maxfaces" - the maximum number of faces per cell
  ! "/face/cell1" and "/face/cell2" - the arrays the face edge data
  module subroutine read_topology(par_env, case_name, topo)

    use constants, only: geoext, adiosconfig
    use io, only: initialise_io, cleanup_io, configure_io, &
                  open_file, close_file, &
                  read_scalar, read_array
    use types, only: io_environment, io_process
    use parallel, only: read_command_line_arguments

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    character(len=:), allocatable :: case_name                              !< The name of the case that is computed
    type(topology), target, intent(inout) :: topo                           !< The topology for which to compute the parition

    ! Local variables
    character(len=:), allocatable :: geo_file    ! Geo file name
    character(len=:), allocatable :: adios2_file ! ADIOS2 config file name

    class(io_environment), allocatable :: io_env
    class(io_process), allocatable :: geo_reader

    integer(ccs_long), dimension(1) :: sel_start
    integer(ccs_long), dimension(1) :: sel_count

    geo_file = case_name//geoext
    adios2_file = case_name//adiosconfig

    call initialise_io(par_env, adios2_file, io_env)
    call configure_io(io_env, "geo_reader", geo_reader)  
  
    call open_file(geo_file, "read", geo_reader)
  
    ! Read attribute "ncel" - the total number of cells
    call read_scalar(geo_reader, "ncel", topo%global_num_cells)
    ! Read attribute "nfac" - the total number of faces
    call read_scalar(geo_reader, "nfac", topo%global_num_faces)
    ! Read attribute "maxfaces" - the maximum number of faces per cell
    call read_scalar(geo_reader, "maxfaces", topo%max_faces)

    allocate(topo%face_edge_end1(topo%global_num_faces))
    allocate(topo%face_edge_end2(topo%global_num_faces))

    sel_start(1) = 0 ! Global index to start reading from
    sel_count(1) = topo%global_num_faces ! How many elements to read in total

    ! Read arrays face/cell1 and face/cell2
    call read_array(geo_reader, "/face/cell1", sel_start, sel_count, topo%face_edge_end1)
    call read_array(geo_reader, "/face/cell2", sel_start, sel_count, topo%face_edge_end2)

    ! Close the file and ADIOS2 engine
    call close_file(geo_reader)

    ! Finalise the ADIOS2 IO environment
    call cleanup_io(io_env)

  end subroutine read_topology

  module subroutine compute_connectivity(par_env, topo)

    use mpi

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(topology), target, intent(inout) :: topo                           !< The topology for which to compute the parition

    ! Local variables
    integer(ccs_int) :: i, j, k, m, n
    integer(ccs_int) :: ierr
    integer(ccs_int) :: num_local_entries ! The number of local P2P entries
    integer(ccs_int) :: num_entries       ! The total number of P2P entries
    integer(ccs_int), dimension(:), allocatable :: tmp_int1d   ! Temporary 1D integer array
    integer(ccs_int), dimension(:,:), allocatable :: tmp_int2d ! Temporary 2D integer array
    integer(ccs_int), dimension(:), allocatable :: p2p_row_idx !< P2P matrix row index
    integer(ccs_int), dimension(:), allocatable :: p2p_col_idx !< P2P matrix col index
    integer(ccs_int), dimension(:), allocatable :: p2p_value   !< P2P matrix values

    allocate(tmp_int1d(par_env%num_procs))
    allocate(tmp_int2d(par_env%num_procs,2))

    tmp_int1d = 0

    associate(num_procs => par_env%num_procs, &
              irank => par_env%proc_id)

      ! All ranks loop over the total number of processes from 0 to num_procs - 1
      do k = 0, num_procs - 1

        tmp_int2d = 0
        
        ! Loop over rank's local vertices
        do i = 1, topo%vtxdist(irank + 2) - topo%vtxdist(irank + 1)

          ! Look up global partition for current vertex "i"
          m = topo%global_partition(i + topo%vtxdist(irank + 1) - 1)

          ! If global partition for vertex "i" equals "k"
          if (m == k) then

            ! Loop over the number of vertices adjacent to vertex "i"
            do j = 1, topo%xadj(i + 1) - topo%xadj(i)

              ! Look up the global partition of the adjacent vertex
              n = topo%global_partition(topo%adjncy(topo%xadj(i) + j - 1))

              ! Is the adjacent vertex on the same partition or not?
              if (m /= n) then
                ! Adjacent vertex on different partition
                tmp_int2d(n, 2) = tmp_int2d(n, 2) + 1
              else
                ! Adjacent vertex on same partition
                tmp_int2d(n, 1) = tmp_int2d(n, 1) + 1
              endif

            enddo

          endif

        enddo

        select type(par_env)
        type is(parallel_environment_mpi)
          call MPI_AllReduce(MPI_IN_PLACE, tmp_int2d, num_procs * 2, MPI_INTEGER, MPI_SUM, par_env%comm, ierr)
        end select

        if (k == irank) then
          tmp_int1d = tmp_int2d( : , 2)
          tmp_int1d(irank + 1) = tmp_int2d(k + 1, 1) 
        endif

      enddo

      ! Count local non-zero entries
      num_local_entries = count(tmp_int1d /= 0)
    
      select type(par_env)
      type is(parallel_environment_mpi)

        ! Compute total number of non-zero entries
        call MPI_Allreduce(num_local_entries, num_entries, 1, MPI_INTEGER, MPI_SUM, par_env%comm, ierr)

      end select

      call dprint("Number of local P2P entries: "//str(num_local_entries))
      call dprint("Total number of P2P entries: "//str(num_entries))

      ! Allocate arrays to store P2P connectivity (as CSR matrix)
      allocate(p2p_value(num_entries))
      allocate(p2p_col_idx(num_entries))
      allocate(p2p_row_idx(num_procs+1))
      
      ! Initialise arrays to 0
      p2p_row_idx = 0
      p2p_col_idx = 0
      p2p_value = 0         

      ! Compute row indices
      p2p_row_idx(irank + 2) = num_local_entries

      select type(par_env)
      type is(parallel_environment_mpi)

        call MPI_AllReduce(MPI_IN_PLACE, p2p_row_idx, num_procs + 1, MPI_INTEGER, MPI_SUM, par_env%comm, ierr)

      end select

      p2p_row_idx(1) = 1
      do j = 2, num_procs + 1
        p2p_row_idx(j) = p2p_row_idx(j) + p2p_row_idx(j-1) 
      enddo
      
      ! Compute column indices and values
      i = p2p_row_idx(irank + 1)

      do j = 1,num_procs
        if (tmp_int1d(j) /= 0) then
          p2p_col_idx(i) = j
          p2p_value(i) = tmp_int1d(j)
          i = i + 1
        endif
      enddo

      select type(par_env)
      type is(parallel_environment_mpi)
        call MPI_AllReduce(MPI_IN_PLACE, p2p_col_idx, num_entries, MPI_INTEGER, MPI_SUM, par_env%comm, ierr)
        call MPI_AllReduce(MPI_IN_PLACE, p2p_value, num_entries, MPI_INTEGER, MPI_SUM, par_env%comm, ierr)
      end select

      do i = 1, num_entries
        call dprint("P2P matrix value "//str(i)//": "//str(int(p2p_value(i))))
        call dprint("P2P matrix col index "//str(i)//": "//str(int(p2p_col_idx(i))))
      end do

      do i = 1, num_procs + 1
        call dprint("P2P matrix row index "//str(i)//": "//str(int(p2p_row_idx(i))))
      end do
      
    end associate

    deallocate(p2p_col_idx)
    deallocate(p2p_row_idx)
    deallocate(p2p_value)
    deallocate(tmp_int1d)
    deallocate(tmp_int2d)

  end subroutine compute_connectivity

end submodule