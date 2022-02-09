module io

  use types, only: io_environment, io_process
  use kinds, only: accs_int, accs_real
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
    module procedure read_scalar_integer
    module procedure read_scalar_real
  end interface

  interface read_array
    module procedure read_array_integer1D
    module procedure read_array_integer2D
    module procedure read_array_real1D
    module procedure read_array_real2D
  end interface

  interface 

  !> @brief Initialise the IO environment
  module subroutine initialise_io(par_env, config_file, io_env)
    class(parallel_environment), intent(in) :: par_env
    character (len=*), optional, intent(in) :: config_file
    class(io_environment), allocatable, intent(out) :: io_env
  end subroutine

  !> @brief Clean up the IO environment
  module subroutine cleanup_io(io_env)
    class(io_environment), intent(inout) :: io_env
  end subroutine

  !> @brief Configure the IO process
  module subroutine configure_io(io_env, process_name, io_proc)
    class(io_environment), intent(in) :: io_env
    character (len=*), intent(in) :: process_name
    class(io_process), allocatable, intent(out) :: io_proc
  end subroutine

  !> brief Open file
  !
  !> param[in] filename : name of file to open
  !> param[in] mode : choose whether to read/ write or append
  !> param[inout] io_handle : object that include IO environment handles
  module subroutine open_file(filename, mode, io_proc)
    character (len=*), intent(in) :: filename
    character (len=*), intent(in) :: mode
    class(io_process), intent(inout) :: io_proc
  end subroutine

  !> brief Close file/engine with ADIOS2
  !
  !> param[in] io_process : ADIOS2 IO process
  module subroutine close_file(io_proc)
    class(io_process), intent(inout) :: io_proc
  end subroutine
  

  module subroutine read_scalar_integer(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: attr_name
    integer(accs_int), intent(out) :: attr
  end subroutine

  module subroutine read_scalar_real(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: attr_name
    real(accs_real), intent(out) :: attr
  end subroutine

  module subroutine read_array_integer1D(io_proc, var_name, start, count, var)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: var_name
    integer(kind=8), dimension(1), intent(in) :: start
    integer(kind=8), dimension(1), intent(in) :: count
    integer, dimension(:), intent(inout) :: var
  end subroutine

  module subroutine read_array_integer2D(io_proc, var_name, start, count, var)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: var_name
    integer(kind=8), dimension(2), intent(in) :: start
    integer(kind=8), dimension(2), intent(in) :: count
    integer, dimension(:,:), intent(inout) :: var
  end subroutine

  module subroutine read_array_real1D(io_proc, var_name, start, count, var)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: var_name
    integer(kind=8), dimension(1), intent(in) :: start
    integer(kind=8), dimension(1), intent(in) :: count
    real, dimension(:), intent(inout) :: var
  end subroutine

  module subroutine read_array_real2D(io_proc, var_name, start, count, var)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: var_name
    integer(kind=8), dimension(2), intent(in) :: start
    integer(kind=8), dimension(2), intent(in) :: count
    real, dimension(:,:), intent(inout) :: var
  end subroutine

  end interface

end module io