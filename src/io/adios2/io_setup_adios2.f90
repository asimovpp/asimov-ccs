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
  use timestepping, only: get_timestep, get_current_time, timestepping_is_active

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

    real(ccs_real) :: sim_time
    real(ccs_real) :: delta_t

    integer(ccs_int) :: ierr
    type(adios2_attribute) :: time_attr
    type(adios2_attribute) :: dt_attr

    allocate (adios2_io_process :: io_proc)

    select type (io_env)
    type is (adios2_env)

      select type (io_proc)
      type is (adios2_io_process)

        call adios2_declare_io(io_proc%io_task, io_env%adios, process_name, ierr)
        if(timestepping_is_active()) then
          delta_t = get_timestep()
          call adios2_define_attribute(dt_attr, io_proc%io_task, "dt", delta_t, ierr)
        end if
        call get_current_time(sim_time)
        call adios2_define_attribute(time_attr, io_proc%io_task, "simulation time", sim_time, ierr)

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

    character(len=:), allocatable :: engine_type
    character(len=:), allocatable :: file_type
    
    integer(ccs_int) :: ierr

    file_type = ""

    select type (io_proc)
    type is (adios2_io_process)

      if(mode == "write") then

        ! query the engine type - defined in the ADIOS2 XML configuration file
        call adios2_io_engine_type(engine_type, io_proc%io_task, ierr)

        ! Support for HDF5, BP4 and BP5
        if (engine_type == "HDF5") then
          file_type = ".h5"
        else if (engine_type == "BP4" .or. engine_type == "BP5") then
          file_type = ".bp"
        else
          call error_abort("Unknown ADIOS2 engine type: " // trim(engine_type))
        end if

        ! Append the correct file extension
        call adios2_open(io_proc%engine, io_proc%io_task, filename // file_type, get_mode(mode), ierr)

      else if(mode == "read" .or. mode == "append") then

        call adios2_open(io_proc%engine, io_proc%io_task, filename, get_mode(mode), ierr)

      end if

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
