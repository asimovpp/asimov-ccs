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
  use types, only: ccs_mesh, io_environment, io_process
  use mesh_utils, only: read_mesh
  use partitioning, only: compute_partitioner_input, &
                          partition_kway, compute_connectivity

  implicit none

  type(ccs_mesh) :: mesh

  class(*), pointer :: config_file_pointer  !< Pointer to CCS config file
  character(len=error_length) :: error

  character(len=:), allocatable :: case_name   !< Case name
  character(len=:), allocatable :: ccs_config_file  !< Config file for CCS
  character(len=:), allocatable :: geo_file
  character(len=:), allocatable :: adios2_file

  class(parallel_environment), allocatable :: par_env

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

  call read_mesh(par_env, case_name, mesh)
  call partition_kway(par_env, mesh)
  call compute_connectivity(par_env, mesh)

  call timer(end_time)

  if (par_env%proc_id == 0) then
    print *, "Elapsed time: ", end_time - start_time
  end if

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
