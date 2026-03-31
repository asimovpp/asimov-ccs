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
  use logging, only: log_unit_out

  implicit none

contains

  !v Creates a new parallel environment by splitting the existing one, splitting
  !  based on provided MPI constants or a provided colouring
  module subroutine create_new_par_env(parent_par_env, split, use_mpi_splitting, par_env)
    class(parallel_environment), intent(in) :: parent_par_env         !< The parent parallel environment
    integer, intent(in) :: split                                      !< The value indicating which type of split is being performed, or the user provided colour
    logical, intent(in) :: use_mpi_splitting                          !< Flag indicating whether to use mpi_comm_split_type
    class(parallel_environment), allocatable, intent(out) :: par_env  !< The resulting parallel environment

    integer :: newcomm
    integer :: colour
    integer(ccs_err) :: ierr

    allocate(parallel_environment_mpi :: par_env)

    select type (parent_par_env)
    type is (parallel_environment_mpi)
      call set_colour_from_split(parent_par_env, split, use_mpi_splitting, colour)
      if (use_mpi_splitting) then
        call mpi_comm_split_type(parent_par_env%comm, colour, 0, MPI_INFO_NULL, newcomm, ierr) 
      else 
        call mpi_comm_split(parent_par_env%comm, colour, 0, newcomm, ierr) 
      end if
      call error_handling(ierr, "mpi", parent_par_env)

      select type (par_env)
      type is (parallel_environment_mpi)
        call create_parallel_environment_from_comm(newcomm, par_env)
      class default
        call error_abort("Unsupported parallel environment")
      end select

    class default
      call error_abort("Unsupported parallel environment")
    end select
  end subroutine create_new_par_env
	
  !> Creates a parallel environment based on the provided communicator
  subroutine create_parallel_environment_from_comm(comm, par_env)
    integer, intent(in) :: comm                                          !< The communicator with which to make the parallel environment
    type(parallel_environment_mpi), intent(inout) :: par_env   !< The resulting parallel environment

    par_env%comm = comm
    call set_mpi_parameters(par_env)
  end subroutine create_parallel_environment_from_comm

  !> Sets mpi parameters inside a parallel environment
  module subroutine set_mpi_parameters(par_env)
    class(parallel_environment), intent(inout) :: par_env !< The parallel environment being updated

    integer(ccs_err) :: ierr

    select type (par_env)
    type is (parallel_environment_mpi)
      if (is_valid(par_env)) then
        call mpi_comm_rank(par_env%comm, par_env%proc_id, ierr)
        call error_handling(ierr, "mpi", par_env)

        call mpi_comm_size(par_env%comm, par_env%num_procs, ierr)
        call error_handling(ierr, "mpi", par_env)

        call par_env%set_rop()
        par_env%root=0
      else
        par_env%proc_id = -1
        par_env%num_procs = 0
        par_env%root=-1
      end if
    class default
      call error_abort("Unsupported parallel environment")
    end select
  end subroutine set_mpi_parameters

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
    integer(mpi_address_kind) :: byte_size, allocate_byte_size

    disp_unit = c_sizeof(dummy_int)

    select type (shared_env)
    type is (parallel_environment_mpi)
      allocate_byte_size = shared_alloc_size([length], disp_unit, shared_env, .true.)
      call mpi_win_allocate_shared(allocate_byte_size, disp_unit, MPI_INFO_NULL, shared_env%comm, c_array_ptr, window, ierr)
      call error_handling(ierr, "mpi", shared_env)

      call mpi_win_shared_query(window, 0, byte_size, disp_unit, c_array_ptr, ierr)
      call error_handling(ierr, "mpi", shared_env)

      call c_f_pointer(c_array_ptr, array, shape=[length])

      call mpi_barrier(shared_env%comm, ierr)
      call error_handling(ierr, "mpi", shared_env)

    class default
      call error_abort("Unsupported parallel environment")

    end select

  end subroutine

  module subroutine create_shared_array_long_1D(shared_env, length, array, window)

    use iso_c_binding

    class(parallel_environment), intent(in) :: shared_env
    integer(ccs_int), intent(in) :: length
    integer(ccs_long), pointer, dimension(:), intent(out) :: array
    integer, intent(out) :: window
    type(c_ptr) :: c_array_ptr
    integer(ccs_long) :: dummy_long = 1_ccs_long
    integer(ccs_err) :: ierr
    integer :: disp_unit
    integer(mpi_address_kind) :: byte_size, allocate_byte_size

    disp_unit = c_sizeof(dummy_long)

    select type (shared_env)
    type is (parallel_environment_mpi)
      allocate_byte_size = shared_alloc_size([length], disp_unit, shared_env, .true.)
      call mpi_win_allocate_shared(allocate_byte_size, disp_unit, MPI_INFO_NULL, shared_env%comm, c_array_ptr, window, ierr)
      call error_handling(ierr, "mpi", shared_env)

      call mpi_win_shared_query(window, 0, byte_size, disp_unit, c_array_ptr, ierr)
      call error_handling(ierr, "mpi", shared_env)

      call c_f_pointer(c_array_ptr, array, shape=[length])

      call mpi_barrier(shared_env%comm, ierr)
      call error_handling(ierr, "mpi", shared_env)

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
    integer(mpi_address_kind) :: byte_size, allocate_byte_size

    disp_unit = c_sizeof(dummy_int)

    select type (shared_env)
    type is (parallel_environment_mpi)
      allocate_byte_size = shared_alloc_size(length, disp_unit, shared_env, .true.)
      call mpi_win_allocate_shared(allocate_byte_size, disp_unit, MPI_INFO_NULL, shared_env%comm, c_array_ptr, window, ierr)
      call error_handling(ierr, "mpi", shared_env)

      call mpi_win_shared_query(window, 0, byte_size, disp_unit, c_array_ptr, ierr)
      call error_handling(ierr, "mpi", shared_env)

      call c_f_pointer(c_array_ptr, array, shape=length)

      call mpi_barrier(shared_env%comm, ierr)
      call error_handling(ierr, "mpi", shared_env)

    class default
      call error_abort("Unsupported parallel environment")

    end select

  end subroutine

  module subroutine create_shared_array_real_1D(shared_env, length, array, window)

    use iso_c_binding

    class(parallel_environment), intent(in) :: shared_env
    integer(ccs_int), intent(in) :: length
    real(ccs_real), pointer, dimension(:), intent(out) :: array
    integer, intent(out) :: window

    type(c_ptr) :: c_array_ptr
    real(ccs_real) :: dummy_real = 1.0_ccs_real
    integer(ccs_err) :: ierr
    integer :: disp_unit
    integer(mpi_address_kind) :: byte_size, allocate_byte_size

    disp_unit = c_sizeof(dummy_real)

    select type (shared_env)
    type is (parallel_environment_mpi)
      allocate_byte_size = shared_alloc_size([length], disp_unit, shared_env, .true.)
      call mpi_win_allocate_shared(allocate_byte_size, disp_unit, MPI_INFO_NULL, shared_env%comm, c_array_ptr, window, ierr)
      call error_handling(ierr, "mpi", shared_env)

      call mpi_win_shared_query(window, 0, byte_size, disp_unit, c_array_ptr, ierr)
      call error_handling(ierr, "mpi", shared_env)

      call c_f_pointer(c_array_ptr, array, shape=[length])

      call mpi_barrier(shared_env%comm, ierr)
      call error_handling(ierr, "mpi", shared_env)

    class default
      call error_abort("Unsupported parallel environment")

    end select

  end subroutine

  module subroutine create_shared_array_real_2D(shared_env, length, array, window)

    use iso_c_binding

    class(parallel_environment), intent(in) :: shared_env
    integer(ccs_int), dimension(2), intent(in) :: length
    real(ccs_real), pointer, dimension(:,:), intent(out) :: array
    integer, intent(out) :: window

    type(c_ptr) :: c_array_ptr
    real(ccs_real) :: dummy_real = 1.0_ccs_real
    integer(ccs_err) :: ierr
    integer :: disp_unit
    integer(mpi_address_kind) :: byte_size, allocate_byte_size

    disp_unit = c_sizeof(dummy_real)

    select type (shared_env)
    type is (parallel_environment_mpi)
      allocate_byte_size = shared_alloc_size(length, disp_unit, shared_env, .true.)
      call mpi_win_allocate_shared(allocate_byte_size, disp_unit, MPI_INFO_NULL, shared_env%comm, c_array_ptr, window, ierr)
      call error_handling(ierr, "mpi", shared_env)

      call mpi_win_shared_query(window, 0, byte_size, disp_unit, c_array_ptr, ierr)
      call error_handling(ierr, "mpi", shared_env)

      call c_f_pointer(c_array_ptr, array, shape=length)

      call mpi_barrier(shared_env%comm, ierr)
      call error_handling(ierr, "mpi", shared_env)

    class default
      call error_abort("Unsupported parallel environment")

    end select

  end subroutine

  module subroutine create_shared_array_real_3D(shared_env, length, array, window)

    use iso_c_binding

    class(parallel_environment), intent(in) :: shared_env
    integer(ccs_int), dimension(3), intent(in) :: length
    real(ccs_real), pointer, dimension(:,:,:), intent(out) :: array
    integer, intent(out) :: window

    type(c_ptr) :: c_array_ptr
    real(ccs_real) :: dummy_real = 1.0_ccs_real
    integer(ccs_err) :: ierr
    integer :: disp_unit
    integer(mpi_address_kind) :: byte_size, allocate_byte_size

    disp_unit = c_sizeof(dummy_real)

    select type (shared_env)
    type is (parallel_environment_mpi)
      allocate_byte_size = shared_alloc_size([length], disp_unit, shared_env, .true.)
      call mpi_win_allocate_shared(allocate_byte_size, disp_unit, MPI_INFO_NULL, shared_env%comm, c_array_ptr, window, ierr)
      call error_handling(ierr, "mpi", shared_env)

      call mpi_win_shared_query(window, 0, byte_size, disp_unit, c_array_ptr, ierr)
      call error_handling(ierr, "mpi", shared_env)

      call c_f_pointer(c_array_ptr, array, shape=length)

      call mpi_barrier(shared_env%comm, ierr)
      call error_handling(ierr, "mpi", shared_env)

    class default
      call error_abort("Unsupported parallel environment")

    end select

  end subroutine

  integer(MPI_ADDRESS_KIND) function shared_alloc_size(shape, elem_bytes, shared_env, root_alloc) result(nbytes)

    use iso_fortran_env, only: int64
    
    integer(ccs_int), dimension(:), intent(in) :: shape      !< What is the shape we are allocating
    integer(ccs_int), intent(in) :: elem_bytes               !< Size of elements (in bytes)
    type(parallel_environment_mpi), intent(in) :: shared_env !< The shared memory environment
    logical, intent(in) :: root_alloc                        !< Option to set allocation pattern

    integer(MPI_ADDRESS_KIND) :: mybytes
    integer(MPI_ADDRESS_KIND) :: delta

    if (any(shape < 0)) then
       call error_abort("Allocated shape must be positive")
    end if
    if (elem_bytes < 0) then
       call error_abort("Per element size [bytes] must be positive")
    end if
    
    ! Determine how many bytes to allocate
    nbytes = product(int(shape, int64)) * int(elem_bytes, int64)

    ! Determine how to distribute the allocation
    if (root_alloc) then
       ! Root rank allocates
       if (.not. is_root(shared_env)) then
          nbytes = 0
       end if
    else
       ! Allocate across shared memory region
       mybytes = nbytes / shared_env%num_procs

       ! Distribute bytes to lower ranks
       delta = nbytes - shared_env%num_procs * mybytes
       if (shared_env%proc_id < delta) then
          mybytes = mybytes + 1
       end if

       nbytes = mybytes
    end if
    
    if (nbytes < 0) then
       call error_abort("Byte computation overflowed!")
    end if
    
  end function shared_alloc_size

  module subroutine destroy_shared_array_int_1D(shared_env, array, window)
    class(parallel_environment), intent(in) :: shared_env
    integer(ccs_int), pointer, dimension(:), intent(inout) :: array
    integer, intent(inout) :: window
    integer(ccs_err) :: ierr

    ! Keeping shared_env as argument to more clearly decrate this function as MPI shared memory related
    associate(foo => shared_env)
    end associate

    call mpi_win_free(window, ierr)
    call error_handling(ierr, "mpi", shared_env)

    nullify(array)

  end subroutine

  module subroutine destroy_shared_array_long_1D(shared_env, array, window)
    class(parallel_environment), intent(in) :: shared_env
    integer(ccs_long), pointer, dimension(:), intent(inout) :: array
    integer, intent(inout) :: window
    integer(ccs_err) :: ierr

    ! Keeping shared_env as argument to more clearly decrate this function as MPI shared memory related
    associate(foo => shared_env)
    end associate

    call mpi_win_free(window, ierr)
    call error_handling(ierr, "mpi", shared_env)

    nullify(array)

  end subroutine

  module subroutine destroy_shared_array_real_1D(shared_env, array, window)
    class(parallel_environment), intent(in) :: shared_env
    real(ccs_real), pointer, dimension(:), intent(inout) :: array
    integer, intent(inout) :: window

    integer(ccs_err) :: ierr

    ! Keeping shared_env as argument to more clearly decrate this function as MPI shared memory related
    associate(foo => shared_env)
    end associate

    call mpi_win_free(window, ierr)
    call error_handling(ierr, "mpi", shared_env)

    nullify(array)

  end subroutine


  module subroutine destroy_shared_array_int_2D(shared_env, array, window)
    class(parallel_environment), intent(in) :: shared_env
    integer(ccs_int), pointer, dimension(:,:), intent(inout) :: array
    integer, intent(inout) :: window
    integer(ccs_err) :: ierr

    ! Keeping shared_env as argument to more clearly decrate this function as MPI shared memory related
    associate(foo => shared_env)
    end associate

    call mpi_win_free(window, ierr)
    call error_handling(ierr, "mpi", shared_env)

    nullify(array)

  end subroutine

  module subroutine destroy_shared_array_real_2D(shared_env, array, window)
    class(parallel_environment), intent(in) :: shared_env
    real(ccs_real), pointer, dimension(:,:), intent(inout) :: array
    integer, intent(inout) :: window
    integer(ccs_err) :: ierr

    ! Keeping shared_env as argument to more clearly decrate this function as MPI shared memory related
    associate(foo => shared_env)
    end associate

    call mpi_win_free(window, ierr)
    call error_handling(ierr, "mpi", shared_env)

    nullify(array)

  end subroutine

  !> Synchronise the parallel environment
  module subroutine sync(par_env)

    class(parallel_environment), intent(in) :: par_env

    integer(ccs_err) :: ierr ! Error code

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

        if (is_root(par_env)) then
          write(log_unit_out,*) new_line('a') // "Usage: ./ccs_app [OPTIONS]" // new_line('a')
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
              if (is_root(par_env)) then
                call print_help()
              end if
              call cleanup_parallel_environment(par_env)
              stop 0
            case default
              if (is_root(par_env)) then
                write(log_unit_out,*) "Argument ", trim(arg), " not supported by ASiMoV-CCS."
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

    write(log_unit_out,*) "========================================="
    write(log_unit_out,*) "ASiMoV-CCS command line OPTIONS          "
    write(log_unit_out,*) "========================================="
    write(log_unit_out,*) "--ccs_help:               This help menu"
    write(log_unit_out,*) "--ccs_m <value>:          Problem size"
    write(log_unit_out,*) "--ccs_case <string>:      Test case name" // new_line('a')
    write(log_unit_out,*) "--ccs_in <string>:        Path to input directory" // new_line('a')

  end subroutine

  !> Query whether a STOP file exists and broadcast result to all processes
  module function query_stop_run(par_env) result(stop_run)

    class(parallel_environment), intent(in) :: par_env !< parallel_environment_mpi
    logical :: stop_run
    integer(ccs_int) :: ierr

    stop_run = .false.

    select type (par_env)
    type is (parallel_environment_mpi)
      if(is_root(par_env)) then 
        inquire(file="STOP", EXIST=stop_run)
        if (stop_run) then
          write(log_unit_out,*) "STOP file found"
        end if
      end if

      call MPI_Bcast(stop_run, 1, MPI_LOGICAL, par_env%root, par_env%comm, ierr)
      call error_handling(ierr, "mpi", par_env)

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

  !> Check whether current process is root process in communicator
  module function is_valid(par_env) result(isvalid)
    class(parallel_environment), intent(in) :: par_env !< parallel environment
    logical :: isvalid

    isvalid = .false.

    select type (par_env)
    type is (parallel_environment_mpi)
      if (par_env%comm == MPI_COMM_NULL) then
        isvalid = .false.
      else if (par_env%comm /= MPI_COMM_NULL) then
        isvalid = .true. 
      else 
        call error_abort("communicator not initialised")
      end if

    class default
      call error_abort("Unsupported parallel environment")

    end select

  end function is_valid

  !> Sets the colour for splitting the parallel environment based on the split value provided
  module subroutine set_colour_from_split(par_env, split_type, use_mpi_splitting, colour)
    use constants

    class(parallel_environment), intent(in) :: par_env    !< The parallel environment
    integer, intent(in) :: split_type                     !< Split value provided
    logical, intent(in) :: use_mpi_splitting              !< Flag indicating whether to use mpi_comm_split_type
    integer, intent(out) :: colour                        !< The resulting colour

    select type (par_env)
    type is (parallel_environment_mpi)
      if (use_mpi_splitting) then
        if (split_type == ccs_split_type_shared) then
          colour = MPI_COMM_TYPE_SHARED
        end if
      else
        if (split_type == ccs_split_undefined) then 
          colour = MPI_UNDEFINED
        else if (split_type == ccs_split_type_low_high) then
          colour = 0
          if (par_env%proc_id >= par_env%num_procs/2) then
            colour = 1
          end if
        else if (split_type >= 0) then
          colour = split_type
        end if
      end if
    class default
      call error_abort("Unsupported parallel environment")
    end select
  end subroutine set_colour_from_split

  !> Creates communicator of roots of specified shared environments
  module subroutine create_shared_roots_comm(par_env, shared_env, roots_env)
    use constants
    class(parallel_environment), intent(in) :: par_env                     !< The parent parallel environment of the shared_envs
    class(parallel_environment), intent(in) :: shared_env                  !< The shared environments whose roots we want in the root environment
    class(parallel_environment), allocatable, intent(out) :: roots_env   !< The resulting root environment

    integer :: colour
    logical :: use_mpi_splitting
    
    if (is_root(shared_env)) then
      colour = 1
    else 
      colour = ccs_split_undefined
    end if

    use_mpi_splitting = .false.

    call create_new_par_env(par_env, colour, use_mpi_splitting, roots_env)

  end subroutine create_shared_roots_comm

end submodule parallel_utils_mpi
