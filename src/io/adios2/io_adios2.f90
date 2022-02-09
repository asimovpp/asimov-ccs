submodule (io) io_adios2

  use adios2
  use adios2_types, only: adios2_env, adios2_io_process
  use parallel_types_mpi, only: parallel_environment_mpi

  implicit none

  contains

  module subroutine read_scalar_integer(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: attr_name
    integer(accs_int), intent(out) :: attr

    type(adios2_attribute) :: adios2_attr

    integer(accs_int) :: ierr

    select type(io_proc)
      type is(adios2_io_process)

        call adios2_inquire_attribute(adios2_attr, io_proc%io_task, attr_name, ierr)
        call adios2_attribute_data(attr, adios2_attr, ierr)

      class default
        print*,"Unknown IO process handler type"

      end select

  end subroutine

  module subroutine read_scalar_real(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: attr_name
    real(accs_real), intent(out) :: attr

    type(adios2_attribute) :: adios2_attr

    integer(accs_int) :: ierr

    select type(io_proc)
      type is(adios2_io_process)

        call adios2_inquire_attribute(adios2_attr, io_proc%io_task, attr_name, ierr)
        call adios2_attribute_data(attr, adios2_attr, ierr)

      class default
        print*,"Unknown IO process handler type"

      end select

    end subroutine

    module subroutine read_array_integer1D(io_proc, var_name, start, count, var)
      class(io_process), intent(in) :: io_proc
      character (len=*), intent(in) :: var_name
      integer(kind=8), dimension(1), intent(in) :: start
      integer(kind=8), dimension(1), intent(in) :: count
      integer, dimension(:), intent(inout) :: var

      type(adios2_variable):: adios2_var
      integer(accs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
          call adios2_set_selection(adios2_var, 1, start, count, ierr)
          call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

      class default
        print*,"Unknown IO process handler type"

      end select

    end subroutine

    module subroutine read_array_integer2D(io_proc, var_name, start, count, var)
      class(io_process), intent(in) :: io_proc
      character (len=*), intent(in) :: var_name
      integer(kind=8), dimension(2), intent(in) :: start
      integer(kind=8), dimension(2), intent(in) :: count
      integer, dimension(:,:), intent(inout) :: var

      type(adios2_variable):: adios2_var
      integer(accs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
          call adios2_set_selection(adios2_var, 2, start, count, ierr)
          call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

      class default
        print*,"Unknown IO process handler type"

      end select

    end subroutine

    module subroutine read_array_real1D(io_proc, var_name, start, count, var)
      class(io_process), intent(in) :: io_proc
      character (len=*), intent(in) :: var_name
      integer(kind=8), dimension(1), intent(in) :: start
      integer(kind=8), dimension(1), intent(in) :: count
      real, dimension(:), intent(inout) :: var

      type(adios2_variable):: adios2_var
      integer(accs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
          call adios2_set_selection(adios2_var, 1, start, count, ierr)
          call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

      class default
        print*,"Unknown IO process handler type"

      end select

    end subroutine

    module subroutine read_array_real2D(io_proc, var_name, start, count, var)
      class(io_process), intent(in) :: io_proc
      character (len=*), intent(in) :: var_name
      integer(kind=8), dimension(2), intent(in) :: start
      integer(kind=8), dimension(2), intent(in) :: count
      real, dimension(:,:), intent(inout) :: var

      type(adios2_variable):: adios2_var
      integer(accs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
          call adios2_set_selection(adios2_var, 2, start, count, ierr)
          call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

      class default
        print*,"Unknown IO process handler type"

      end select

    end subroutine


  end submodule