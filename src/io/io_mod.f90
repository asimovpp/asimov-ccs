!v Module file io_mod.f90
!
!  Provides an interface to IO functions.
module io

  use iso_fortran_env, only: int32, int64, real32, real64
  use types, only: io_environment, io_process
  use parallel_types, only: parallel_environment
  use constants, only: ndim, adiosconfig
  use kinds, only: ccs_int, ccs_real, ccs_long

  implicit none

  private

  public :: initialise_io
  public :: cleanup_io
  public :: configure_io
  public :: open_file
  public :: close_file
  public :: get_num_steps
  public :: read_scalar
  public :: read_array
  public :: write_scalar
  public :: write_array

  interface read_scalar
    module procedure read_scalar_int32
    module procedure read_scalar_int64
    module procedure read_scalar_real32
    module procedure read_scalar_real64
  end interface

  interface read_array
    module procedure read_array_int32_1D
    module procedure read_array_int64_1D
    module procedure read_array_int32_2D
    module procedure read_array_int64_2D
    module procedure read_array_real32_1D
    module procedure read_array_real64_1D
    module procedure read_array_real32_2D
    module procedure read_array_real64_2D
  end interface

  interface write_scalar
    module procedure write_scalar_int32
    module procedure write_scalar_int64
    module procedure write_scalar_real32
    module procedure write_scalar_real64
  end interface

  interface write_array
    module procedure write_array_int32_1D
    module procedure write_array_int64_1D
    module procedure write_array_int32_2D
    module procedure write_array_int64_2D
    module procedure write_array_real32_1D
    module procedure write_array_real64_1D
    module procedure write_array_real32_2D
    module procedure write_array_real64_2D
  end interface

  interface

    !> Initialise the IO environment
    module subroutine initialise_io(par_env, config_file, io_env)
      class(parallel_environment), intent(in) :: par_env        !< parallel environment that IO environment will reside on
      character(len=*), optional, intent(in) :: config_file     !< name of the IO configuration file
      class(io_environment), allocatable, intent(out) :: io_env !< IO environment
    end subroutine

    !> Clean up the IO environment
    module subroutine cleanup_io(io_env)
      class(io_environment), intent(inout) :: io_env !< IO environment
    end subroutine

    !> Configure the IO process
    module subroutine configure_io(io_env, process_name, io_proc)
      class(io_environment), intent(in) :: io_env            !< IO environment
      character(len=*), intent(in) :: process_name           !< name of the IO process to be configured
      class(io_process), allocatable, intent(out) :: io_proc !< the configured IO process
    end subroutine

    !> Open file
    module subroutine open_file(filename, mode, io_proc)
      character(len=*), intent(in) :: filename    !< name of file to open
      character(len=*), intent(in) :: mode        !< choose whether to read, write or append
      class(io_process), intent(inout) :: io_proc !< object that include IO environment handles
    end subroutine

    !> Close file
    module subroutine close_file(io_proc)
      class(io_process), intent(inout) :: io_proc !< IO process
    end subroutine

    module subroutine get_num_steps(io_proc, steps)
      class(io_process), intent(inout) :: io_proc
      integer(ccs_long), intent(out) :: steps
    end subroutine

    !> Read a scalar integer from file
    module subroutine read_scalar_int32(io_proc, attr_name, attr)
      class(io_process), intent(in) :: io_proc  !< IO process used for reading
      character(len=*), intent(in) :: attr_name !< Name of scalar integer to read
      integer(int32), intent(out) :: attr       !< Value of scalar integer
    end subroutine

    !> Read a scalar long integer from file
    module subroutine read_scalar_int64(io_proc, attr_name, attr)
      class(io_process), intent(in) :: io_proc  !< IO process used for reading
      character(len=*), intent(in) :: attr_name !< Name of scalar long integer to read
      integer(int64), intent(out) :: attr       !< Value of scalar long integer
    end subroutine

    !> Read a scalar real from file
    module subroutine read_scalar_real32(io_proc, attr_name, attr)
      class(io_process), intent(in) :: io_proc  !< IO process used for reading
      character(len=*), intent(in) :: attr_name !< Name of scalar real to read
      real(real32), intent(out) :: attr         !< Value of scalar real
    end subroutine

    !> Read a scalar double precision real from file
    module subroutine read_scalar_real64(io_proc, attr_name, attr)
      class(io_process), intent(in) :: io_proc  !< IO process used for reading
      character(len=*), intent(in) :: attr_name !< Name of scalar double precision real to read
      real(real64), intent(out) :: attr         !< Value of scalar double precision real
    end subroutine

    !> Read a 1D 32-bit integer array from file
    module subroutine read_array_int32_1D(io_proc, var_name, global_start, count, var, step)
      class(io_process), intent(in) :: io_proc                 !< IO process used for reading
      character(len=*), intent(in) :: var_name                 !< Name of integer array to read
      integer(int64), dimension(1), intent(in) :: global_start !< What global index to start reading from
      integer(int64), dimension(1), intent(in) :: count        !< How many array element to read
      integer(int32), dimension(:), intent(inout) :: var       !< The 1D integer array
      integer(int64), optional, intent(in) :: step             !< The step to read
    end subroutine

    !> Read a 1D 64-bit integer array from file
    module subroutine read_array_int64_1D(io_proc, var_name, global_start, count, var, step)
      class(io_process), intent(in) :: io_proc                 !< IO process used for reading
      character(len=*), intent(in) :: var_name                 !< Name of integer array to read
      integer(int64), dimension(1), intent(in) :: global_start !< What global index to start reading from
      integer(int64), dimension(1), intent(in) :: count        !< How many array element to read
      integer(int64), dimension(:), intent(inout) :: var       !< The 1D integer array
      integer(int64), optional, intent(in) :: step             !< The step to read
    end subroutine

    !> Read a 2D 32-bit integer array from file
    module subroutine read_array_int32_2D(io_proc, var_name, global_start, count, var, step)
      class(io_process), intent(in) :: io_proc                 !< IO process used for reading
      character(len=*), intent(in) :: var_name                 !< Name of integer array to read
      integer(int64), dimension(2), intent(in) :: global_start !< What global index to start reading from
      integer(int64), dimension(2), intent(in) :: count        !< How many array element to read
      integer(int32), dimension(:, :), intent(inout) :: var    !< The 2D integer array
      integer(int64), optional, intent(in) :: step             !< The step to read
    end subroutine

    !> Read a 2D 64-bit integer array from file
    module subroutine read_array_int64_2D(io_proc, var_name, global_start, count, var, step)
      class(io_process), intent(in) :: io_proc                 !< IO process used for reading
      character(len=*), intent(in) :: var_name                 !< Name of integer array to read
      integer(int64), dimension(2), intent(in) :: global_start !< What global index to start reading from
      integer(int64), dimension(2), intent(in) :: count        !< How many array element to read
      integer(int64), dimension(:, :), intent(inout) :: var    !< The 2D integer array
      integer(int64), optional, intent(in) :: step             !< The step to read
    end subroutine

    !> Read a 1D 32-bit real array from file
    module subroutine read_array_real32_1D(io_proc, var_name, global_start, count, var, step)
      class(io_process), intent(in) :: io_proc                 !< IO process used for reading
      character(len=*), intent(in) :: var_name                 !< Name of real array to read
      integer(int64), dimension(1), intent(in) :: global_start !< What global index to start reading from
      integer(int64), dimension(1), intent(in) :: count        !< How many array element to read
      real(real32), dimension(:), intent(inout) :: var         !< The 1D real array
      integer(int64), optional, intent(in) :: step             !< The step to read
    end subroutine

    !> Read a 1D 64-bit real array from file
    module subroutine read_array_real64_1D(io_proc, var_name, global_start, count, var, step)
      class(io_process), intent(in) :: io_proc                 !< IO process used for reading
      character(len=*), intent(in) :: var_name                 !< Name of real array to read
      integer(int64), dimension(1), intent(in) :: global_start !< What global index to start reading from
      integer(int64), dimension(1), intent(in) :: count        !< How many array element to read
      real(real64), dimension(:), intent(inout) :: var         !< The 1D real array
      integer(int64), optional, intent(in) :: step             !< The step to read
    end subroutine

    !> Read a 2D 32-bit real array from file
    module subroutine read_array_real32_2D(io_proc, var_name, global_start, count, var, step)
      class(io_process), intent(in) :: io_proc                 !< IO process used for reading
      character(len=*), intent(in) :: var_name                 !< Name of real array to read
      integer(int64), dimension(2), intent(in) :: global_start !< What global index to start reading from
      integer(int64), dimension(2), intent(in) :: count        !< How many array element to read
      real(real32), dimension(:, :), intent(inout) :: var      !< The 2D real array
      integer(int64), optional, intent(in) :: step             !< The step to read
    end subroutine

    !> Read a 2D 64-bit real array from file
    module subroutine read_array_real64_2D(io_proc, var_name, global_start, count, var, step)
      class(io_process), intent(in) :: io_proc                 !< IO process used for reading
      character(len=*), intent(in) :: var_name                 !< Name of real array to read
      integer(int64), dimension(2), intent(in) :: global_start !< What global index to start reading from
      integer(int64), dimension(2), intent(in) :: count        !< How many array element to read
      real(real64), dimension(:, :), intent(inout) :: var      !< The 2D real array
      integer(int64), optional, intent(in) :: step             !< The step to read
    end subroutine

    !>  Write a scalar integer to file
    module subroutine write_scalar_int32(io_proc, attr_name, attr)
      class(io_process), intent(in) :: io_proc    !< IO process used for reading
      character(len=*), intent(in) :: attr_name  !< Name of scalar integer to read
      integer(int32), intent(in) :: attr          !< Value of scalar integer
    end subroutine

    !>  Write a scalar long integer to file
    module subroutine write_scalar_int64(io_proc, attr_name, attr)
      class(io_process), intent(in) :: io_proc   !< IO process used for reading
      character(len=*), intent(in) :: attr_name !< Name of scalar integer to read
      integer(int64), intent(in) :: attr         !< Value of scalar integer
    end subroutine

    !>  Write a scalar real to file
    module subroutine write_scalar_real32(io_proc, attr_name, attr)
      class(io_process), intent(in) :: io_proc   !< IO process used for reading
      character(len=*), intent(in) :: attr_name !< Name of scalar real to read
      real(real32), intent(in) :: attr           !< Value of scalar real
    end subroutine

    !>  Write a scalar double precision real to file
    module subroutine write_scalar_real64(io_proc, attr_name, attr)
      class(io_process), intent(in) :: io_proc   !< IO process used for reading
      character(len=*), intent(in) :: attr_name !< Name of scalar double precision real to read
      real(real64), intent(in) :: attr           !< Value of scalar double precision real
    end subroutine

    !>  Write a 1D 32-bit integer array to file
    module subroutine write_array_int32_1D(io_proc, var_name, global_shape, global_start, count, var)
      class(io_process), intent(in) :: io_proc                 !< IO process used for reading
      character(len=*), intent(in) :: var_name                !< Name of integer array to read
      integer(int64), dimension(1), intent(in) :: global_shape !< Global shape of array
      integer(int64), dimension(1), intent(in) :: global_start !< The global index to start reading from
      integer(int64), dimension(1), intent(in) :: count        !< How many array elements to read
      integer(int32), dimension(:), intent(in) :: var          !< The 1D integer array
    end subroutine

    !>  Write a 1D 64-bit integer array to file
    module subroutine write_array_int64_1D(io_proc, var_name, global_shape, global_start, count, var)
      class(io_process), intent(in) :: io_proc                 !< IO process used for reading
      character(len=*), intent(in) :: var_name                !< Name of integer array to read
      integer(int64), dimension(1), intent(in) :: global_shape !< Global shape of array
      integer(int64), dimension(1), intent(in) :: global_start !< The global index to start reading from
      integer(int64), dimension(1), intent(in) :: count        !< How many array elements to read
      integer(int64), dimension(:), intent(in) :: var          !< The 1D integer array
    end subroutine

    !>  Write a 2D 32-bit integer array to file
    module subroutine write_array_int32_2D(io_proc, var_name, global_shape, global_start, count, var)
      class(io_process), intent(in) :: io_proc                 !< IO process used for reading
      character(len=*), intent(in) :: var_name                !< Name of integer array to read
      integer(int64), dimension(2), intent(in) :: global_shape !< Global shape of array
      integer(int64), dimension(2), intent(in) :: global_start !< The global index to start reading from
      integer(int64), dimension(2), intent(in) :: count        !< How many array elements to read
      integer(int32), dimension(:, :), intent(in) :: var        !< The 2D integer array
    end subroutine

    !>  Write a 2D 64-bit integer array to file
    module subroutine write_array_int64_2D(io_proc, var_name, global_shape, global_start, count, var)
      class(io_process), intent(in) :: io_proc                 !< IO process used for reading
      character(len=*), intent(in) :: var_name                !< Name of integer array to read
      integer(int64), dimension(2), intent(in) :: global_shape !< Global shape of array
      integer(int64), dimension(2), intent(in) :: global_start !< The global index to start reading from
      integer(int64), dimension(2), intent(in) :: count        !< How many array elements to read
      integer(int64), dimension(:, :), intent(in) :: var        !< The 2D integer array
    end subroutine

    !>  Write a 1D 32-bit real array to file
    module subroutine write_array_real32_1D(io_proc, var_name, global_shape, global_start, count, var)
      class(io_process), intent(in) :: io_proc                 !< IO process used for reading
      character(len=*), intent(in) :: var_name                !< Name of real array to read
      integer(int64), dimension(1), intent(in) :: global_shape !< Global shape of array
      integer(int64), dimension(1), intent(in) :: global_start !< The global index to start reading from
      integer(int64), dimension(1), intent(in) :: count        !< How many array elements to read
      real(real32), dimension(:), intent(in) :: var            !< The 1D real array
    end subroutine

    !>  Write a 1D 64-bit real array to file
    module subroutine write_array_real64_1D(io_proc, var_name, global_shape, global_start, count, var)
      class(io_process), intent(in) :: io_proc                 !< IO process used for reading
      character(len=*), intent(in) :: var_name                !< Name of real array to read
      integer(int64), dimension(1), intent(in) :: global_shape !< Global shape of array
      integer(int64), dimension(1), intent(in) :: global_start !< The global index to start reading from
      integer(int64), dimension(1), intent(in) :: count        !< How many array elements to read
      real(real64), dimension(:), intent(inout) :: var         !< The 1D real array
    end subroutine

    !>  Write a 2D 32-bit real array to file
    module subroutine write_array_real32_2D(io_proc, var_name, global_shape, global_start, count, var)
      class(io_process), intent(in) :: io_proc                 !< IO process used for reading
      character(len=*), intent(in) :: var_name                !< Name of real array to read
      integer(int64), dimension(2), intent(in) :: global_shape !< Global shape of array
      integer(int64), dimension(2), intent(in) :: global_start !< The global index to start reading from
      integer(int64), dimension(2), intent(in) :: count        !< How many array elements to read
      real(real32), dimension(:, :), intent(in) :: var          !< The 2D real array
    end subroutine

    !>  Write a 2D 64-bit real array to file
    module subroutine write_array_real64_2D(io_proc, var_name, global_shape, global_start, count, var)
      class(io_process), intent(in) :: io_proc                 !< IO process used for reading
      character(len=*), intent(in) :: var_name                !< Name of real array to read
      integer(int64), dimension(2), intent(in) :: global_shape !< Global shape of arra
      integer(int64), dimension(2), intent(in) :: global_start !< The global index to start reading from
      integer(int64), dimension(2), intent(in) :: count        !< How many array elements to read
      real(real64), dimension(:, :), intent(in) :: var          !< The 2D real array
    end subroutine

  end interface

end module io
