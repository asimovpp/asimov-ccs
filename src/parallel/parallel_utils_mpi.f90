!v Submodule file parallel_utils_mpi.smod
!
!  Implementation (using MPI) of the parallel utilities
!
!  @build mpi

submodule(parallel) parallel_utils_mpi
#include "ccs_macros.inc"

  use utils, only: exit_print
  use mpi
  use parallel_types_mpi, only: parallel_environment_mpi

  implicit none

contains

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
  module subroutine read_command_line_arguments(par_env, cps, case_name)

    class(parallel_environment), intent(in) :: par_env !< parallel_environment_mpi
    integer(ccs_int), optional, intent(inout) :: cps   !< number of cells per side
    character(len=:), optional, allocatable, intent(out) :: case_name

    character(len=1024) :: arg ! argument string
    integer(ccs_int) :: nargs  ! number of arguments
    integer(ccs_int) :: arg_len

    arg_len = 0

    select type (par_env)
    type is (parallel_environment_mpi)

      if (command_argument_count() == 0) then

        if (par_env%proc_id == par_env%root) then
          print *, "========================================="
          print *, "ASiMoV-CCS needs command line options:"
          print *, "--ccs_help:               The help menu."
          print *, "--ccs_m <value>:          Problem size."
          print *, "--ccs_case <string>:      Test case name."
          print *, "========================================="
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
              call get_command_argument(nargs + 1, length = arg_len, value = arg)
              if (present(case_name)) then
                allocate (character(len = arg_len) :: case_name)
                case_name = trim(arg)
              end if
            case ('--ccs_help')
              if (par_env%proc_id == par_env%root) then
                print *, "========================================="
                print *, "ASiMoV-CCS command line options:         "
                print *, "========================================="
                print *, "--ccs_help:               This help menu."
                print *, "--ccs_m <value>:          Problem size."
                print *, "--ccs_case <string>:      Test case name."
              end if
              call cleanup_parallel_environment(par_env)
              stop 0
            case default
              if (par_env % proc_id == par_env % root) then
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

end submodule parallel_utils_mpi
