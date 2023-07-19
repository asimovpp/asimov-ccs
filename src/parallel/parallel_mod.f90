!v Module file parallel.mod
!
!  Module that defines the parallel interace for ASiMoV-CCS

module parallel

  use parallel_types
  use kinds, only: ccs_int

  implicit none

  private

  public :: initialise_parallel_environment
  public :: create_new_par_env
  public :: cleanup_parallel_environment
  public :: sync
  public :: read_command_line_arguments
  public :: timer
  public :: allreduce
  public :: error_handling !TODO: consider if this should be public (used with "raw" MPI calls in some places)
  public :: query_stop_run
  public :: create_shared_array
  public :: is_root

  interface

    !> Create the parallel environment
    module subroutine initialise_parallel_environment(par_env)
      class(parallel_environment), allocatable, intent(out) :: par_env
    end subroutine

    !> Creates a new parallel environment by splitting the existing one, splitting
    !  based on provided MPI constants or a provided colouring
    module subroutine create_new_par_env(parent_par_env, split, use_mpi_splitting, par_env)
      class(parallel_environment), intent(in) :: parent_par_env         !< The parent parallel environment
      integer, intent(in) :: split                                      !< The value indicating which type of split is being performed, or the user provided colour
      logical, intent(in) :: use_mpi_splitting                          !< Flag indicating whether to use mpi_comm_split_type
      class(parallel_environment), allocatable, intent(out) :: par_env  !< The resulting parallel environment
    end subroutine


    !> Cleanup the parallel environment
    module subroutine cleanup_parallel_environment(par_env)
      class(parallel_environment), intent(in) :: par_env
    end subroutine

    !> Synchronise the parallel environment
    module subroutine sync(par_env)
      class(parallel_environment), intent(in) :: par_env
    end subroutine

    !> Read command line arguments and their values
    module subroutine read_command_line_arguments(par_env, cps, case_name, in_dir)
      class(parallel_environment), intent(in) :: par_env
      integer(ccs_int), optional, intent(inout) :: cps
      character(len=:), optional, allocatable, intent(out) :: case_name
      character(len=:), optional, allocatable, intent(out) :: in_dir
    end subroutine read_command_line_arguments

    !> Timer for parallel environment
    module subroutine timer(tick)
      double precision, intent(out) :: tick
    end subroutine

    !> Global reduction of integer scalars
    module subroutine allreduce_scalar(input_value, rop, par_env, result_value)
      class(*), intent(in) :: input_value
      class(reduction_operator), intent(in) :: rop
      class(parallel_environment), intent(in) :: par_env
      class(*), intent(inout) :: result_value
    end subroutine

    !> Error handling for parallel environment
    module subroutine error_handling(error_code, error_category, par_env)
      integer, intent(in) :: error_code
      character(len=*), intent(in) :: error_category
      class(parallel_environment), intent(in) :: par_env
    end subroutine

    !> Query whether a STOP file exists
    module function query_stop_run(par_env) result(stop_run)
      class(parallel_environment), intent(in) :: par_env !< parallel_environment_mpi
      logical :: stop_run
    end function

    !> Create an integer 1D MPI shared memory array
    module subroutine create_shared_array_int_1D(shared_env, length, array, window)
      class(parallel_environment), intent(in) :: shared_env
      integer(ccs_int), intent(in) :: length
      integer(ccs_int), pointer, dimension(:), intent(out) :: array
      integer, intent(out) :: window
    end subroutine

    !> Create an integer 2D MPI shared memory array
    module subroutine create_shared_array_int_2D(shared_env, length, array, window)
      class(parallel_environment), intent(in) :: shared_env
      integer(ccs_int), dimension(2), intent(in) :: length
      integer(ccs_int), pointer, dimension(:,:), intent(out) :: array
      integer, intent(out) :: window
    end subroutine

    !> Test whether current rank is root of communicator
    module function is_root(par_env) result(isroot)
      class(parallel_environment), intent(in) :: par_env
      logical :: isroot
    end function

  end interface

  interface create_shared_array
    module procedure create_shared_array_int_1D
    module procedure create_shared_array_int_2D
  end interface

  interface allreduce
    module procedure allreduce_scalar
  end interface allreduce

end module parallel
