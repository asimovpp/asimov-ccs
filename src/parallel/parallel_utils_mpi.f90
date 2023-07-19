!v Submodule file parallel_utils_mpi.smod
!
!  Implementation (using MPI) of the parallel utilities
!
!  @build mpi
!  @dont_fail_linter

submodule(parallel) parallel_utils_mpi
#include "ccs_macros.inc"

  use utils, only: exit_print
  use mpi
  use parallel_types_mpi, only: parallel_environment_mpi
  use kinds, only: ccs_err

  implicit none

contains


  module subroutine create_shared_array_int_1D(shared_env, length, array, window)

    use iso_c_binding

    class(parallel_environment), intent(in) :: shared_env
    integer(ccs_int), intent(in) :: length
    integer(ccs_int), pointer, dimension(:), intent(out) :: array
    integer, intent(out) :: window
    type(c_ptr) :: c_array_ptr
    integer(ccs_int) :: dummy_int = 1_ccs_int
    integer(ccs_err) :: ierr
    integer :: disp_unit
    integer(mpi_address_kind) :: base_ptr, byte_size, allocate_byte_size

    disp_unit = c_sizeof(dummy_int)
    byte_size = length * disp_unit

    if (isroot(shared_env)) then
      allocate_byte_size = byte_size
    else
      allocate_byte_size = 0
    end if

    select type (shared_env)
    type is (parallel_environment_mpi)
      call mpi_win_allocate_shared(allocate_byte_size, disp_unit, MPI_INFO_NULL, shared_env%comm, c_array_ptr, window, ierr)

      call c_f_pointer(c_array_ptr, array, shape=[length])

      call mpi_barrier(shared_env%comm, ierr)

      call mpi_win_shared_query(window, 0, byte_size, disp_unit, base_ptr, ierr)

    class default
      call error_abort("Unsupported parallel environment")

    end select

  end subroutine

  module subroutine create_shared_array_int_2D(shared_env, length, array, window)

    use iso_c_binding

    class(parallel_environment), intent(in) :: shared_env
    integer(ccs_int), dimension(2), intent(in) :: length
    integer(ccs_int), pointer, dimension(:,:), intent(out) :: array
    integer, intent(out) :: window
    type(c_ptr) :: c_array_ptr
    integer(ccs_int) :: dummy_int = 1_ccs_int
    integer(ccs_err) :: ierr
    integer :: disp_unit
    integer(mpi_address_kind) :: base_ptr, byte_size, allocate_byte_size

    disp_unit = c_sizeof(dummy_int)
    byte_size = length(1) * length(2) * disp_unit

    if (isroot(shared_env)) then
      allocate_byte_size = byte_size
    else
      allocate_byte_size = 0
    end if

    select type (shared_env)
    type is (parallel_environment_mpi)
      call mpi_win_allocate_shared(allocate_byte_size, disp_unit, MPI_INFO_NULL, shared_env%comm, c_array_ptr, window, ierr)

      call c_f_pointer(c_array_ptr, array, shape=[length])

      call mpi_barrier(shared_env%comm, ierr)

      call mpi_win_shared_query(window, 0, byte_size, disp_unit, base_ptr, ierr)

    class default
      call error_abort("Unsupported parallel environment")

    end select

  end subroutine

  !> Synchronise the parallel environment
  module subroutine sync(par_env)

    class(parallel_environment), intent(in) :: par_env

    integer :: ierr ! Error code

    select type (par_env)
    type is (parallel_environment_mpi)

      call mpi_barrier(par_env%comm, ierr)
      call error_handling(ierr, "mpi", par_env)

    class default
      call error_abort("Unsupported parallel environment")

    end select

  end subroutine

  !> read command line arguments and their values
  module subroutine read_command_line_arguments(par_env, cps, case_name, in_dir)

    class(parallel_environment), intent(in) :: par_env !< parallel_environment_mpi
    integer(ccs_int), optional, intent(inout) :: cps   !< number of cells per side
    character(len=:), optional, allocatable, intent(out) :: case_name
    character(len=:), optional, allocatable, intent(out) :: in_dir

    character(len=1024) :: arg ! argument string
    integer(ccs_int) :: nargs  ! number of arguments
    integer(ccs_int) :: arg_len

    arg_len = 0

    select type (par_env)
    type is (parallel_environment_mpi)

      if (command_argument_count() == 0) then

        if (par_env%proc_id == par_env%root) then
          print *, new_line('a') // "Usage: ./ccs_app [OPTIONS]" // new_line('a')
          call print_help()
          call cleanup_parallel_environment(par_env)
          stop 0
        end if

      else

        do nargs = 1, command_argument_count()
          call get_command_argument(nargs, arg)

          if (arg(1:6) == '--ccs_') then
            select case (arg)
            case ('--ccs_m') ! problems size
              call get_command_argument(nargs + 1, arg)
              read (arg, '(I5)') cps
            case ('--ccs_case') ! case name
              call get_command_argument(nargs + 1, length=arg_len, value=arg)
              if (present(case_name)) then
                allocate (character(len=arg_len) :: case_name)
                case_name = trim(arg)
              end if
            case ('--ccs_in') ! input directory
              call get_command_argument(nargs + 1, length=arg_len, value=arg)
              if (present(in_dir)) then
                allocate (character(len=arg_len) :: in_dir)
                in_dir = trim(arg)
              end if
            case ('--ccs_help')
              if (par_env%proc_id == par_env%root) then
                call print_help()
              end if
              call cleanup_parallel_environment(par_env)
              stop 0
            case default
              if (par_env%proc_id == par_env%root) then
                print *, "Argument ", trim(arg), " not supported by ASiMoV-CCS."
              end if
              call cleanup_parallel_environment(par_env)
              stop 1
            end select
          end if
        end do

      end if

    class default
      call error_abort("Unsupported parallel environment")
    end select

  end subroutine read_command_line_arguments

  !> Timer for MPI parallel environment
  module subroutine timer(tick)

    double precision, intent(out) :: tick !< variable that is assigned the current time

    tick = mpi_wtime()

  end subroutine

  subroutine print_help()

    print *, "========================================="
    print *, "ASiMoV-CCS command line OPTIONS          "
    print *, "========================================="
    print *, "--ccs_help:               This help menu"
    print *, "--ccs_m <value>:          Problem size"
    print *, "--ccs_case <string>:      Test case name" // new_line('a')
    print *, "--ccs_in <string>:        Path to input directory" // new_line('a')

  end subroutine

    !> Query whether a STOP file exists and broadcast result to all processes
  module function query_stop_run(par_env) result(stop_run)

    class(parallel_environment), intent(in) :: par_env !< parallel_environment_mpi
    logical :: stop_run
    integer(ccs_int) :: ierr

    stop_run = .false.

    select type (par_env)
    type is (parallel_environment_mpi)
      if(par_env%proc_id == par_env%root) then 
        inquire(file="STOP", EXIST=stop_run)
      end if

      call MPI_Bcast(stop_run, 1, MPI_LOGICAL, par_env%root, par_env%comm, ierr)

    class default
      call error_abort("Unsupported parallel environment")
    end select

  end function

  !> Check whether current process is root process in communicator
  module function is_root(par_env) result(isroot)
    class(parallel_environment), intent(in) :: par_env !< parallel environment
    logical :: isroot

    isroot = .false.

    select type (par_env)
    type is (parallel_environment_mpi)
      if (par_env%proc_id == par_env%root) then
        isroot = .true.
      end if

    class default
      call error_abort("Unsupported parallel environment")

    end select

  end function is_root

end submodule parallel_utils_mpi
