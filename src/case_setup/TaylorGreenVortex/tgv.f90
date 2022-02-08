
program tgv

  use, intrinsic :: iso_fortran_env, only:  output_unit

  use yaml, only: parse, error_length
  use read_config, only: get_case_name
  use kinds, only : accs_int, accs_real
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment
  use parallel_types, only: parallel_environment
  use types, only: io_environment, &
                   io_process
  use io, only: initialise_io, cleanup_io, configure_io, &
                open_file, close_file, &
                read_attribute

  implicit none

  class(*), pointer :: config_file
  character(len=error_length) :: error

  ! Case title
  character(len=:), allocatable :: case_name
  ! Geo file name
  character(len=:), allocatable :: geo_file

  class(parallel_environment), allocatable :: par_env
  class(io_environment), allocatable :: io_env
  class(io_process), allocatable :: geo_reader

!  integer(accs_int) :: irank, isize
  integer(accs_int) :: max_faces
  integer(accs_int) :: num_faces
  integer(accs_int) :: num_cells

  ! Launch MPI
  call initialise_parallel_environment(par_env) 

  call initialise_io(par_env, "adios2-config.xml", io_env)

  ! Read case name from configuration file
  call read_configuration()

  geo_file = case_name//".geo"
  
  call configure_io(io_env, "test_reader", geo_reader)

  call open_file(geo_file, "read", geo_reader)

  call read_attribute(geo_reader, "ncel", num_cells)
  call read_attribute(geo_reader, "nfac", num_faces)
  call read_attribute(geo_reader, "maxfaces", max_faces)

  print*, "Max number of faces: ", max_faces
  print*, "Total number of faces: ", num_faces
  print*, "Total number of cells: ", num_cells

  call close_file(geo_reader)

  call cleanup_io(io_env)

  call cleanup_parallel_environment(par_env)

  contains

  subroutine read_configuration()

    config_file => parse("./tgv_config.yaml", error=error)
    if (error/='') then
      print*,trim(error)
      stop 1
    endif
    
    ! Get title
    call get_case_name(config_file, case_name)

  end subroutine

end program tgv
