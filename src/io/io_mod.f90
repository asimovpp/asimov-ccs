module io

  use types, only: io_environment, io_process
  use kinds, only: accs_int
  use parallel_types, only: parallel_environment

  implicit none

  private

  public :: initialise_io
  public :: configure_io
  public :: open_file
  public :: close_file

  interface

  !> @brief Initialise the IO environment
  module subroutine initialise_io(par_env, config_file, io_env)
    class(parallel_environment), intent(in) :: par_env
    character (len=*), allocatable, optional, intent(in) :: config_file
    class(io_environment), intent(inout) :: io_env
  end subroutine

  !> @brief Configure the IO process
  module subroutine configure_io(process_name, io_env, io)
    character (len=*), intent(in) :: process_name
    class(io_environment), intent(inout) :: io_env
    class(io_process), intent(inout) :: io
  end subroutine

  !> brief Open file
  !
  !> param[in] filename : name of file to open
  !> param[in] mode : choose whether to read/ write or append
  !> param[inout] io_handle : object that include IO environment handles
  module subroutine open_file(filename, mode, io)
    character (len=*), intent(in) :: filename
    integer(accs_int), intent(in) :: mode
    class(io_process), intent(inout) :: io
  end subroutine

  module subroutine close_file()
  end subroutine
  
  end interface

end module io