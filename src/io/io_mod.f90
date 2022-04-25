!>  Module file io_mod.f90
!
!>  Provides an interface to IO functions.
module io

  use iso_fortran_env, only: int32, int64, real32, real64
  use types, only: io_environment, io_process
  use parallel_types, only: parallel_environment

  implicit none

  private

  public :: initialise_io
  public :: cleanup_io
  public :: configure_io
  public :: open_file
  public :: close_file
  public :: read_scalar
  public :: read_array

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

  interface 

  !>  Initialise the IO environment
  !
  !> param[in]  par_env     : parallel environment that IO environment
  !!                          will reside on
  !> param[in]  config_file : name of the IO configuration file
  !> param[out] io_env      : IO environment
  module subroutine initialise_io(par_env, config_file, io_env)
    class(parallel_environment), intent(in) :: par_env
    character (len=*), optional, intent(in) :: config_file
    class(io_environment), allocatable, intent(out) :: io_env
  end subroutine

  !>  Clean up the IO environment
  !
  !> param[inout] io_env : IO environment
  module subroutine cleanup_io(io_env)
    class(io_environment), intent(inout) :: io_env
  end subroutine

  !>  Configure the IO process
  !
  !> param[in]  io_env       : IO environment
  !> param[in]  process_name : name of the IO process to be configured
  !> param[out] io_proc      : the configured IO process
  module subroutine configure_io(io_env, process_name, io_proc)
    class(io_environment), intent(in) :: io_env
    character (len=*), intent(in) :: process_name
    class(io_process), allocatable, intent(out) :: io_proc
  end subroutine

  !>  Open file
  !
  !> param[in] filename   : name of file to open
  !> param[in] mode       : choose whether to read, write or append
  !> param[inout] io_proc : object that include IO environment handles
  module subroutine open_file(filename, mode, io_proc)
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: mode
    class(io_process), intent(inout) :: io_proc
  end subroutine

  !>  Close file
  !
  !> param[in] io_proc : IO process
  module subroutine close_file(io_proc)
    class(io_process), intent(inout) :: io_proc
  end subroutine

  !>  Read a scalar integer from file
  !
  !> param[in]  io_proc   : IO process used for reading
  !> param[in]  attr_name : Name of scalar integer to read
  !> param[out] attr      : Value of scalar integer
  module subroutine read_scalar_int32(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: attr_name
    integer(int32), intent(out) :: attr
  end subroutine

    !>  Read a scalar long integer from file
  !
  !> param[in]  io_proc   : IO process used for reading
  !> param[in]  attr_name : Name of scalar long integer to read
  !> param[out] attr      : Value of scalar long integer
  module subroutine read_scalar_int64(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: attr_name
    integer(int64), intent(out) :: attr
  end subroutine

  !>  Read a scalar real from file
  !
  !> param[in]  io_proc   : IO process used for reading
  !> param[in]  attr_name : Name of scalar real to read
  !> param[out] attr      : Value of scalar real
  module subroutine read_scalar_real32(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: attr_name
    real(real32), intent(out) :: attr
  end subroutine

  !>  Read a scalar double precision real from file
  !
  !> param[in]  io_proc   : IO process used for reading
  !> param[in]  attr_name : Name of scalar double precision real to read
  !> param[out] attr      : Value of scalar double precision real
  module subroutine read_scalar_real64(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: attr_name
    real(real64), intent(out) :: attr
  end subroutine

  !>  Read a 1D 32-bit integer array from file
  !
  !> param[in]    io_proc  : IO process used for reading
  !> param[in]    var_name : Name of integer array to read
  !> param[in]    global_start : What global index to start reading from
  !> param[in]    count    : How many array element to read
  !> param[input] var      : The 1D integer array
  module subroutine read_array_int32_1D(io_proc, var_name, global_start, count, var)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: var_name
    integer(int64), dimension(1), intent(in) :: global_start
    integer(int64), dimension(1), intent(in) :: count
    integer(int32), dimension(:), intent(inout) :: var
  end subroutine
  !>  Read a 1D 64-bit integer array from file
  !
  !> param[in]    io_proc  : IO process used for reading
  !> param[in]    var_name : Name of integer array to read
  !> param[in]    global_start : What global index to start reading from
  !> param[in]    count    : How many array element to read
  !> param[input] var      : The 1D integer array
    module subroutine read_array_int64_1D(io_proc, var_name, global_start, count, var)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: var_name
    integer(int64), dimension(1), intent(in) :: global_start
    integer(int64), dimension(1), intent(in) :: count
    integer(int64), dimension(:), intent(inout) :: var
  end subroutine

  !>  Read a 2D 32-bit integer array from file
  !
  !> param[in]    io_proc  : IO process used for reading
  !> param[in]    var_name : Name of integer array to read
  !> param[in]    global_start : What global index to start reading from
  !> param[in]    count    : How many array element to read
  !> param[input] var      : The 2D integer array
  module subroutine read_array_int32_2D(io_proc, var_name, global_start, count, var)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: var_name
    integer(int64), dimension(2), intent(in) :: global_start
    integer(int64), dimension(2), intent(in) :: count
    integer(int32), dimension(:,:), intent(inout) :: var
  end subroutine

  !>  Read a 2D 64-bit integer array from file
  !
  !> param[in]    io_proc  : IO process used for reading
  !> param[in]    var_name : Name of integer array to read
  !> param[in]    global_start : What global index to start reading from
  !> param[in]    count    : How many array element to read
  !> param[input] var      : The 2D integer array
    module subroutine read_array_int64_2D(io_proc, var_name, global_start, count, var)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: var_name
    integer(int64), dimension(2), intent(in) :: global_start
    integer(int64), dimension(2), intent(in) :: count
    integer(int64), dimension(:,:), intent(inout) :: var
  end subroutine

  !>  Read a 1D 32-bit real array from file
  !
  !> param[in]    io_proc  : IO process used for reading
  !> param[in]    var_name : Name of real array to read
  !> param[in]    global_start : What global index to start reading from
  !> param[in]    count    : How many array element to read
  !> param[input] var      : The 1D real array
  module subroutine read_array_real32_1D(io_proc, var_name, global_start, count, var)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: var_name
    integer(int64), dimension(1), intent(in) :: global_start
    integer(int64), dimension(1), intent(in) :: count
    real(real32), dimension(:), intent(inout) :: var
  end subroutine

  !>  Read a 1D 64-bit real array from file
  !
  !> param[in]    io_proc  : IO process used for reading
  !> param[in]    var_name : Name of real array to read
  !> param[in]    global_start : What global index to start reading from
  !> param[in]    count    : How many array element to read
  !> param[input] var      : The 1D real array
  module subroutine read_array_real64_1D(io_proc, var_name, global_start, count, var)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: var_name
    integer(int64), dimension(1), intent(in) :: global_start
    integer(int64), dimension(1), intent(in) :: count
    real(real64), dimension(:), intent(inout) :: var
  end subroutine

  !>  Read a 2D 32-bit real array from file
  !
  !> param[in]    io_proc  : IO process used for reading
  !> param[in]    var_name : Name of real array to read
  !> param[in]    global_start : What global index to start reading from
  !> param[in]    count    : How many array element to read
  !> param[input] var      : The 2D real array
  module subroutine read_array_real32_2D(io_proc, var_name, global_start, count, var)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: var_name
    integer(int64), dimension(2), intent(in) :: global_start
    integer(int64), dimension(2), intent(in) :: count
    real(real32), dimension(:,:), intent(inout) :: var
  end subroutine

  !>  Read a 2D 64-bit real array from file
  !
  !> param[in]    io_proc  : IO process used for reading
  !> param[in]    var_name : Name of real array to read
  !> param[in]    global_start : What global index to start reading from
  !> param[in]    count    : How many array element to read
  !> param[input] var      : The 2D real array
  module subroutine read_array_real64_2D(io_proc, var_name, global_start, count, var)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: var_name
    integer(int64), dimension(2), intent(in) :: global_start
    integer(int64), dimension(2), intent(in) :: count
    real(real64), dimension(:,:), intent(inout) :: var
  end subroutine

  end interface

end module io
