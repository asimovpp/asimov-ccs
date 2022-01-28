submodule (io) io_setup_adios2

  use adios2
  use adios2_types, only: adios2_env, adios2_io_process
  use parallel_types_mpi, only: parallel_environment_mpi

  implicit none

  contains

  module subroutine initialise_io(par_env, config_file, io_env)
    class(parallel_environment), intent(in) :: par_env
    character (len=*), allocatable, optional, intent(in) :: config_file
    class(io_environment), intent(inout) :: io_env

    integer(accs_int) :: ierr

    select type(io_env)
    type is(adios2_env)

      if(present(config_file)) then
        call adios2_init(io_env%adios, trim(config_file), par_env%comm, adios2_debug_mode_on, ierr)
      else
        print*,"ADIOS2 requires a config file!"
      endif

    class default
      print*, "Unknown IO environment"

    end select

  end subroutine

  module subroutine configure_io(process_name, io_env, io)
    class(io_process), intent(inout) :: io
    class(io_environment), intent(inout) :: io_env
    character (len=*), intent(in) :: process_name

    integer(accs_int) :: ierr

    select type(io_process)
    type is(adios2_io_process)

      call adios2_declare_io(io%io, io_env%adios, process_name, ierr)
      
    class default

      print*,"Unknow IO process handler type"

    end select

  end subroutine

  !> brief Open file with ADIOS2
  !
  !> param[in] filename : name of file to open
  !> param[in] mode : choose whether to read/ write or append
  !!                  valid options are: 
  !!                  adios2_mode_write, adios2_mode_append, adios2_mode_read
  !> param[inout] io_handle : object that include ADIOS2 handler information
  module subroutine open_file(filename, mode, io)
    character (len=*), intent(in) :: filename
    integer(accs_int) :: mode
    class(io_process), intent(inout) :: io

    integer(accs_int) :: ierr

    select type(io_process)
    type is(adios2_io_process)

      call adios2_open(io%engine, io%io, filename, mode, ierr)
      
    class default

      print*,"Unknow IO process handler type"

    end select

  end subroutine

end submodule