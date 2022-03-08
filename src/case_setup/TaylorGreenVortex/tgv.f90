!> @brief Program file for TaylorGreenVortex case
program tgv

  use, intrinsic :: iso_fortran_env, only:  output_unit

  use yaml, only: parse, error_length
  use read_config, only: get_case_name
  use kinds, only : accs_int, accs_real
  use constants, only: geoext, adiosconfig, ccsconfig, ndim
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, &
                      read_command_line_arguments, &
                      timer
  use parallel_types, only: parallel_environment
  use types, only: io_environment, &
                   io_process
  use io, only: initialise_io, cleanup_io, configure_io, &
                open_file, close_file, &
                read_scalar, read_array

  implicit none

  class(*), pointer :: config_file_pointer  !> Pointer to CCS config file
  character(len=error_length) :: error

  character(len=:), allocatable :: case_name   !> Case name
  character(len=:), allocatable :: ccs_config_file  !> Config file for CCS
  character(len=:), allocatable :: geo_file    !> Geo file name
  character(len=:), allocatable :: adios2_file !> ADIOS2 config file name

  class(parallel_environment), allocatable :: par_env
  class(io_environment), allocatable :: io_env
  class(io_process), allocatable :: geo_reader

  real(accs_real), dimension(:,:), allocatable :: xyz_coords
  integer(accs_int), dimension(:), allocatable :: vtxdist
  integer(kind=8), dimension(2) :: xyz_sel_start
  integer(kind=8), dimension(2) :: xyz_sel_count

  integer(accs_int) :: irank !> MPI rank ID
  integer(accs_int) :: isize !> Size of MPI world
  integer(accs_int) :: i, j, k  
  integer(accs_int) :: local_idx_start
  integer(accs_int) :: local_idx_end
  integer(accs_int) :: max_faces !> Maximum number of faces per cell
  integer(accs_int) :: num_faces !> Total number of faces
  integer(accs_int) :: num_cells !> Total number of cells

  double precision :: start_time, end_time

  ! Launch MPI
  call initialise_parallel_environment(par_env) 

  irank = par_env%proc_id
  isize = par_env%num_procs

  call read_command_line_arguments(par_env, case_name=case_name)

  ccs_config_file = case_name//ccsconfig

  ! Read case name from configuration file
  call read_configuration()

  geo_file = case_name//geoext
  adios2_file = case_name//adiosconfig

  call timer(start_time)
  
  call initialise_io(par_env, adios2_file, io_env)
  call configure_io(io_env, "geo_reader", geo_reader)  

  call open_file(geo_file, "read", geo_reader)

  ! Read attributes "ncel", nfac" and "maxfaces"
  call read_scalar(geo_reader, "ncel", num_cells)
  call read_scalar(geo_reader, "nfac", num_faces)
  call read_scalar(geo_reader, "maxfaces", max_faces)

  ! Array to store cell range assigned to each process      
  allocate(vtxdist(isize+1))

  ! Set for element to 1 and last element to world size + 1
  vtxdist(1) = 1
  vtxdist(isize + 1) = num_cells + 1

  ! Divide the total number of cells by the world size to
  ! compute the chunk sizes
  k = int(real(num_cells) / isize)
  j = 1
  do i = 1, isize
     vtxdist(i) = j
     j = j + k
  enddo

  ! First and last cell index assigned to this process
  local_idx_start = vtxdist(irank + 1)
  local_idx_end = vtxdist(irank + 2) - 1

  ! Starting point for reading chunk of data
  xyz_sel_start = (/ 0, vtxdist(irank + 1) - 1 /)
  ! How many data points will be read?
  xyz_sel_count = (/ ndim, vtxdist(irank + 2) - vtxdist(irank + 1)/)

  ! Allocate memory for XYZ coordinates array on each MPI rank
  allocate(xyz_coords(xyz_sel_count(1), xyz_sel_count(2)))

  ! Read XYZ coordinates for variable "/cell/x" 
  call read_array(geo_reader, "/cell/x", xyz_sel_start , xyz_sel_count, xyz_coords)

  ! Close the file and ADIOS2 engine
  call close_file(geo_reader)

  ! Finalise the ADIOS2 IO environment
  call cleanup_io(io_env)

  call timer(end_time)

  if(irank == 0) then
     print*, "Elapsed time: ", end_time - start_time
  end if
  
  ! Deallocate memory for XYZ coordinates array
  deallocate(xyz_coords)

  ! Finalise MPI
  call cleanup_parallel_environment(par_env)

  contains

  ! Read YAML configuration file
  subroutine read_configuration()

    config_file_pointer => parse(ccs_config_file, error=error)
    if (error/='') then
      print*,trim(error)
      stop 1
    endif
    
    ! Get case name
    call get_case_name(config_file_pointer, case_name)

  end subroutine

end program tgv
