!> @brief Submodule file parallel_utils_mpi.smod
!
!> @build mpi
!
!> @details Implementation (using MPI) of the parallel utilities

submodule (parallel) parallel_utils_mpi

  use mpi
  use parallel_types_mpi, only: parallel_environment_mpi

  implicit none

  contains

  !> @brief Synchronise the parallel environment
  !
  !> @param[in] parallel_environment_mpi par_env
  module subroutine sync(par_env)
  
    class(parallel_environment), intent(in) :: par_env
    integer :: ierr !> Error code

    select type (par_env)
      type is (parallel_environment_mpi)
         
        call mpi_barrier(par_env%comm, ierr)
        call error_handling(ierr, "mpi", par_env)

      class default
        write(*,*) "Unsupported parallel environment"
        stop 1

    end select

  end subroutine

  !> @brief read command line arguments and their values
  !
  !> @param[in] par_env : parallel_environment_mpi 
  !> @param[inout] cps: number of cells per side 
  module subroutine read_command_line_arguments(par_env, cps)

    class(parallel_environment), intent(in) :: par_env
    integer(accs_int), optional, intent(inout) :: cps

    character(len=32) :: arg !> argument string
    integer(accs_int) :: nargs !> number of arguments

    select type (par_env)
      type is (parallel_environment_mpi)

        do nargs = 1, command_argument_count()
          call get_command_argument(nargs, arg)
    
          if(arg(1:6) == '--ccs_') then
          select case (arg)
            case ('--ccs_m') !> problems size
              call get_command_argument(nargs+1, arg)
              read(arg, '(I5)') cps
            case ('--ccs_help')
              if(par_env%proc_id == par_env%root) then
                print *, "================================"
                print *, "ASiMoV-CCS command line options:"
                print *, "================================"
                print *, "--ccs_help:         This help menu."
                print *, "--ccs_m <value>:    Problem size."
              end if
              call cleanup_parallel_environment(par_env)
              stop
            case default
              if(par_env%proc_id == par_env%root) then
                print *, "Argument ",trim(arg)," not supported by ASiMoV-CCS!"
              end if
              call cleanup_parallel_environment(par_env)
              stop
            end select
          end if
        end do 

      class default
        write(*,*) "Unsupported parallel environment"
        stop 1
    end select

  end subroutine read_command_line_arguments

  !> @brief Timer for MPI parallel environment
  !
  !> @param[inout] double precision tick - Variable that is assigned the current time
  module subroutine timer(tick)

    double precision, intent(out) :: tick

    tick = mpi_wtime()

  end subroutine

end submodule parallel_utils_mpi
