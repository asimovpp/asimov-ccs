program test_tgv_reader

  use kinds, only: accs_int
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment
  use parallel_types, only: parallel_environment
  use types, only: io_environment, &
                   io_process
  use io, only: initialise_io, cleanup_io, configure_io, &
                open_file, close_file, &
                read_scalar
  
  implicit none

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

!  irank = par_env%proc_id
!  isize = par_env%num_procs

  call configure_io(io_env, "test_reader", geo_reader)

  call open_file("TaylorGreen.geo", "read", geo_reader)

  call read_scalar(geo_reader, "ncel", num_cells)
  call read_scalar(geo_reader, "nfac", num_faces)
  call read_scalar(geo_reader, "maxfaces", max_faces)

  print*, "Max number of faces: ", max_faces
  print*, "Total number of faces: ", num_faces
  print*, "Total number of cells: ", num_cells

  call close_file(geo_reader)

  call cleanup_io(io_env)

  call cleanup_parallel_environment(par_env)

end program test_tgv_reader
