!> @brief Submodule file io_setup_adios2.smod
!
!> @build mpi adios2
!
!> @details Implementation (using MPI and ADIOS@) of parallel 
!!          IO setup functionality
submodule (io) io_setup_adios2

  use adios2
  use adios2_types, only: adios2_env, adios2_io_process
  use parallel_types_mpi, only: parallel_environment_mpi
  use kinds, only: accs_int, accs_real

  implicit none

  contains

  !> @brief Initialise the IO environment
  !
  !> @todo The "mode" is currently hard coded - would be better if this
  !        was a configuration file option?
  !
  !> param[in]  par_env     : parallel environment that IO environment
  !!                          will reside on
  !> param[in]  config_file : name of the ADIOS2 IO configuration file
  !> param[out] io_env      : ADIOS2 IO environment
  module subroutine initialise_io(par_env, config_file, io_env)
    class(parallel_environment), intent(in) :: par_env
    character (len=*), optional, intent(in) :: config_file
    class(io_environment), allocatable, intent(out) :: io_env

    integer(accs_int) :: ierr

    allocate(adios2_env :: io_env)
    
    select type(io_env)
      type is(adios2_env)

        select type(par_env)
          type is(parallel_environment_mpi)

            if(present(config_file)) then
              ! Initialise the ADIOS2 environment
              call adios2_init(io_env%adios, trim(config_file), par_env%comm, adios2_debug_mode_on, ierr)
            else
              print*,"ADIOS2 requires a config file!"
            endif

          class default
            print*, "Unknown parallel environment"
          end select

      class default
        print*, "Unknown IO environment"

    end select

  end subroutine

  !> @brief Clean up the IO environment
  !
  !> param[inout] io_env : ADIOS2 IO environment
  module subroutine cleanup_io(io_env)
    class(io_environment), intent(inout) :: io_env

    integer(accs_int) :: ierr

    select type(io_env)
      type is(adios2_env)

        ! Finalise ADIOS2 environment
        call adios2_finalize(io_env%adios, ierr)

      class default
        print*, "Unknown IO environment"

    end select

  end subroutine
  
  !> @brief Configure the IO process
  !
  !> param[in]  io_env       : ADIOS2 IO environment
  !> param[in]  process_name : name of the IO process to be configured - 
  !!                           must match a name defined in the ADIOS2 
  !!                           configuration XML file
  !> param[out] io_proc      : the configured ADIOS2 IO process
  module subroutine configure_io(io_env, process_name, io_proc)
    class(io_environment), intent(in) :: io_env
    character (len=*), intent(in) :: process_name
    class(io_process), allocatable, intent(out) :: io_proc

    integer(accs_int) :: ierr

    allocate(adios2_io_process :: io_proc)

    select type(io_env)
      type is(adios2_env)

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_declare_io(io_proc%io_task, io_env%adios, process_name, ierr)

        class default

          print*,"Unknown IO process handler type"

      end select

      class default
        print*, "Unknown IO environment"
  
    end select

  end subroutine

  !> @brief Open file with ADIOS2
  !
  !> param[in] filename    : name of file to open
  !> param[in] mode        : choose whether to read/ write or append
  !!                         valid options are: "read", "write", "append"
  !> param[inout] io_proc : object that includes ADIOS2 handler information
  module subroutine open_file(filename, mode, io_proc)
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: mode
    class(io_process), intent(inout) :: io_proc

    integer(accs_int) :: ierr

    select type(io_proc)
      type is(adios2_io_process)

        call adios2_open(io_proc%engine, io_proc%io_task, filename, get_mode(mode), ierr)
      
      class default
        print*,"Unknown IO process handler type"

    end select

  end subroutine

  !> @brief Close file/engine with ADIOS2
  !
  !> param[inout] io_proc : ADIOS2 IO process
  module subroutine close_file(io_proc)
    class(io_process), intent(inout) :: io_proc

    integer(accs_int) :: ierr

    select type(io_proc)
      type is(adios2_io_process)

        call adios2_close(io_proc%engine, ierr)
      
      class default
        print*,"Unknown IO process handler type"

    end select

  end subroutine

  ! @brief Return ADIOS2 values for read, write and append modes
  function get_mode(mode_name) result(mode)
    character (len=*), intent(in) :: mode_name
    integer(accs_int) :: mode

    select case(mode_name)
      case ("read")
        mode = adios2_mode_read
      case ("write")
        mode = adios2_mode_write
      case ("append")
        mode = adios2_mode_append
      case default
        mode = -1
        print*, "Not a valid file mode!"  
      end select

  end function

end submodule