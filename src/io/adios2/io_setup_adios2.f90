!v Submodule file io_setup_adios2.smod
!
!  Implementation (using MPI and ADIOS2) of parallel IO setup functionality
!
!  @build mpi adios2
submodule(io) io_setup_adios2
#include "ccs_macros.inc"

  use utils, only: exit_print
  use adios2
  use adios2_types, only: adios2_env, adios2_io_process
  use parallel_types_mpi, only: parallel_environment_mpi

  implicit none

contains

  !v Initialise the IO environment
  !
  !  @todo The "mode" is currently hard coded - would be better if this
  !        was a configuration file option?
  module subroutine initialise_io(par_env, config_file, io_env)
    class(parallel_environment), intent(in) :: par_env        !< parallel environment that IO environment will reside on
    character(len=*), optional, intent(in) :: config_file     !< name of the ADIOS2 IO configuration file
    class(io_environment), allocatable, intent(out) :: io_env !< ADIOS2 IO environment

    integer(ccs_int) :: ierr

    allocate (adios2_env :: io_env)

    select type (io_env)
    type is (adios2_env)

      select type (par_env)
      type is (parallel_environment_mpi)

        if (present(config_file)) then
          ! Initialise the ADIOS2 environment
          call adios2_init(io_env%adios, trim(config_file), par_env%comm, ierr)
        else
          call error_abort("ADIOS2 requires a config file!")
        end if

      class default
        call error_abort("Unknown parallel environment")
      end select

    class default
      call error_abort("Unknown IO environment")

    end select

  end subroutine

  !> Clean up the IO environment
  module subroutine cleanup_io(io_env)
    class(io_environment), intent(inout) :: io_env !< ADIOS2 IO environment

    integer(ccs_int) :: ierr

    select type (io_env)
    type is (adios2_env)

      ! Finalise ADIOS2 environment
      call adios2_finalize(io_env%adios, ierr)

    class default
      call error_abort("Unknown IO environment")

    end select

  end subroutine

  !> Configure the IO process
  module subroutine configure_io(io_env, process_name, io_proc)
    class(io_environment), intent(in) :: io_env            !< ADIOS2 IO environment
    character(len=*), intent(in) :: process_name           !< name of the IO process to be configured - must match a name
    !< defined in the ADIOS2 configuration XML file
    class(io_process), allocatable, intent(out) :: io_proc !< the configured ADIOS2 IO process

    integer(ccs_int) :: ierr

    allocate (adios2_io_process :: io_proc)

    select type (io_env)
    type is (adios2_env)

      select type (io_proc)
      type is (adios2_io_process)

        call adios2_declare_io(io_proc%io_task, io_env%adios, process_name, ierr)

      class default

        call error_abort("Unknown IO process handler type")

      end select

    class default
      call error_abort("Unknown IO environment")

    end select

  end subroutine

  !> Open file with ADIOS2
  module subroutine open_file(filename, mode, io_proc)
    character(len=*), intent(in) :: filename    !< name of file to open
    character(len=*), intent(in) :: mode        !< choose whether to read/ write or append valid options are:
    !< "read", "write", "append"
    class(io_process), intent(inout) :: io_proc !< object that includes ADIOS2 handler information

    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_open(io_proc%engine, io_proc%io_task, filename, get_mode(mode), ierr)

    class default
      call error_abort("Unknown IO process handler type")

    end select

  end subroutine

  !> Close file/engine with ADIOS2
  module subroutine close_file(io_proc)
    class(io_process), intent(inout) :: io_proc !< ADIOS2 IO process

    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_close(io_proc%engine, ierr)

    class default
      call error_abort("Unknown IO process handler type")

    end select

  end subroutine

  !> Return ADIOS2 values for read, write and append modes
  function get_mode(mode_name) result(mode)
    character(len=*), intent(in) :: mode_name
    integer(ccs_int) :: mode

    select case (mode_name)
    case ("read")
      mode = adios2_mode_read
    case ("write")
      mode = adios2_mode_write
    case ("append")
      mode = adios2_mode_append
    case default
      mode = -1
      call error_abort("Not a valid file mode!")
    end select

  end function

  !> Get the number of steps in the ADIOS2 file
  module subroutine get_num_steps(io_proc, steps)
    class(io_process), intent(inout) :: io_proc
    integer(ccs_long), intent(out) :: steps

    integer(ccs_int) :: ierr

    select type(io_proc)
    type is (adios2_io_process)

      call adios2_steps(steps, io_proc%engine, ierr)

    class default
      call error_abort("Unknown IO process handler type")

    end select

  end subroutine

end submodule
