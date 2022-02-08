submodule (io) io_adios2

  use adios2
  use adios2_types, only: adios2_env, adios2_io_process
  use parallel_types_mpi, only: parallel_environment_mpi

  implicit none

  contains

  module subroutine read_attribute_integer(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: attr_name
    integer(accs_int), intent(out) :: attr

    type(adios2_attribute) :: attr_var

    integer(accs_int) :: ierr

    select type(io_proc)
      type is(adios2_io_process)

        call adios2_inquire_attribute(attr_var, io_proc%io_task, attr_name, ierr)
        call adios2_attribute_data(attr, attr_var, ierr)

      class default
        print*,"Unknown IO process handler type"

      end select

  end subroutine

  module subroutine read_attribute_real(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: attr_name
    real(accs_real), intent(out) :: attr

    type(adios2_attribute) :: attr_var

    integer(accs_int) :: ierr

    select type(io_proc)
      type is(adios2_io_process)

        call adios2_inquire_attribute(attr_var, io_proc%io_task, attr_name, ierr)
        call adios2_attribute_data(attr, attr_var, ierr)

      class default
        print*,"Unknown IO process handler type"

      end select

    end subroutine

  end submodule