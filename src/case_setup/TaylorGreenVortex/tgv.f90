!> Program file for TaylorGreenVortex case
program tgv
#include "ccs_macros.inc"

  use utils, only: exit_print
  use yaml, only: parse, error_length
  use read_config, only: get_case_name
  use constants, only: ccsconfig
  use kinds, only: ccs_int, ccs_real, ccs_long
  use constants, only: geoext, adiosconfig, ccsconfig, ndim
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, &
                      read_command_line_arguments, &
                      timer
  use parallel_types, only: parallel_environment
  use types, only: topology, io_environment, io_process
  use partitioning, only: read_topology, compute_partitioner_input, &
                          partition_kway, compute_connectivity


  implicit none

  type(topology) :: topo

  class(*), pointer :: config_file_pointer  !< Pointer to CCS config file
  character(len=error_length) :: error

  character(len=:), allocatable :: case_name   !< Case name
  character(len=:), allocatable :: ccs_config_file  !< Config file for CCS
  character(len=:), allocatable :: geo_file
  character(len=:), allocatable :: adios2_file

  class(parallel_environment), allocatable :: par_env
  class(io_environment), allocatable :: io_env
  class(io_process), allocatable :: geo_reader

  real(ccs_real), dimension(:, :), allocatable :: xyz_coords
  integer(ccs_long), dimension(2) :: xyz_sel_start
  integer(ccs_long), dimension(2) :: xyz_sel_count

  double precision :: start_time, end_time

  ! Launch MPI
  call initialise_parallel_environment(par_env)

  call read_command_line_arguments(par_env, case_name=case_name)

  ccs_config_file = case_name // ccsconfig

  ! Read case name from configuration file
  call read_configuration()

  geo_file = case_name // geoext
  adios2_file = case_name // adiosconfig

  call timer(start_time)

  call read_topology(par_env, case_name, topo)
  call compute_partitioner_input(par_env, topo)
  call partition_kway(par_env, topo)
  call compute_connectivity(par_env, topo)

  call timer(end_time)

  if(par_env%proc_id == 0) then
    print*, "Elapsed time: ", end_time - start_time
  end if

  ! Starting point for reading chunk of data
  xyz_sel_start = (/0, int(topo%vtxdist(par_env%proc_id + 1)) - 1/)
  ! How many data points will be read?
  xyz_sel_count = (/ndim, int(topo%vtxdist(par_env%proc_id + 2) - topo%vtxdist(par_env%proc_id + 1))/)

  ! Allocate memory for XYZ coordinates array on each MPI rank
  allocate (xyz_coords(xyz_sel_count(1), xyz_sel_count(2)))

  ! Read XYZ coordinates for variable "/cell/x"
  call initialise_io(par_env, adios2_file, io_env)
  call configure_io(io_env, "geo_reader", geo_reader)  
  call open_file(geo_file, "read", geo_reader)
  call read_array(geo_reader, "/cell/x", xyz_sel_start, xyz_sel_count, xyz_coords)
  call close_file(geo_reader) ! Close the file and ADIOS2 engine

  ! Finalise the ADIOS2 IO environment
  call cleanup_io(io_env)

  call timer(end_time)

  if (irank == 0) then
    print *, "Elapsed time: ", end_time - start_time
  end if

  ! Deallocate memory for XYZ coordinates array
  deallocate (xyz_coords)

  ! Finalise MPI
  call cleanup_parallel_environment(par_env)

contains

  ! Read YAML configuration file
  subroutine read_configuration()

    config_file_pointer => parse(ccs_config_file, error=error)
    if (error /= '') then
      call error_abort(trim(error))
    end if

    ! Get case name
    call get_case_name(config_file_pointer, case_name)

  end subroutine

end program tgv
