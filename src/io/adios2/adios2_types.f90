module adios2_types

  use adios2
  use types, only: io_environment, io_process
  use kinds, only: accs_int

  implicit none

  private

  type, public, extends(io_environment) :: adios2_env
    type(adios2_adios):: adios
  end type

  type, public, extends(io_process) :: adios2_io_process
    type(adios2_io):: io_task
    type(adios2_engine):: engine
  end type


end module