!v Submodule file io_adios2.smod
!
! Implementation (using MPI and ADIOS2) of parallel IO functionality
!
! @build mpi adios2
submodule(io) io_adios2
#include "ccs_macros.inc"

  use utils, only: exit_print
  use adios2
  use adios2_types, only: adios2_env, adios2_io_process
  use kinds, only: ccs_long

  implicit none

contains

  !> Read a scalar 32-bit integer from file
  module subroutine read_scalar_int32(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc  !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: attr_name !< Name of scalar integer to read
    integer(int32), intent(out) :: attr       !< Value of scalar integer

    type(adios2_attribute) :: adios2_attr

    integer(int64) :: tmp_attr64
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_attribute(adios2_attr, io_proc%io_task, attr_name, ierr)

      if (adios2_attr%type == adios2_type_integer8) then

        call downcast_warning()
        call adios2_attribute_data(tmp_attr64, adios2_attr, ierr)
        attr = int(tmp_attr64, int32)

      else if (adios2_attr%type == adios2_type_integer4) then

        call adios2_attribute_data(attr, adios2_attr, ierr)

      else
        call error_abort("IO Error: unsuported integer type")
      end if

    class default
      call error_abort("Unknown IO process handler type")

    end select

  end subroutine

  !> Read a scalar 64-bit integer from file
  module subroutine read_scalar_int64(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc  !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: attr_name !< Name of scalar longinteger to read
    integer(int64), intent(out) :: attr       !< Value of scalar long integer

    type(adios2_attribute) :: adios2_attr

    integer(int32) :: tmp_attr32
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_attribute(adios2_attr, io_proc%io_task, attr_name, ierr)

      if (adios2_attr%type == adios2_type_integer4) then

        call adios2_attribute_data(tmp_attr32, adios2_attr, ierr)
        attr = int(tmp_attr32, int64)

      else if (adios2_attr%type == adios2_type_integer8) then

        call adios2_attribute_data(attr, adios2_attr, ierr)

      else
        call error_abort("IO Error: unsuported integer type")
      end if

    class default
      call error_abort("Unknown IO process handler type")

    end select

  end subroutine

  !> Read a scalar 32-bit real from file
  module subroutine read_scalar_real32(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc  !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: attr_name !< Name of scalar real to read
    real(real32), intent(out) :: attr         !< Value of scalar real

    type(adios2_attribute) :: adios2_attr

    real(real64) :: tmp_attr64
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_attribute(adios2_attr, io_proc%io_task, attr_name, ierr)
      if (adios2_attr%type == adios2_type_dp) then

        call downcast_warning()
        call adios2_attribute_data(tmp_attr64, adios2_attr, ierr)
        attr = real(tmp_attr64, real32)

      else if (adios2_attr%type == adios2_type_real) then

        call adios2_attribute_data(attr, adios2_attr, ierr)

      else
        call error_abort("IO Error: unsuported real type")
      end if

    class default
      call error_abort("Unknown IO process handler type")

    end select

  end subroutine

  !> Read a scalar 64-bit real from file
  module subroutine read_scalar_real64(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc  !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: attr_name !< Name of scalar double precision real to read
    real(real64), intent(out) :: attr         !< Value of scalar double precision real

    type(adios2_attribute) :: adios2_attr

    real(real32) :: tmp_attr32
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_attribute(adios2_attr, io_proc%io_task, attr_name, ierr)
      if (adios2_attr%type == adios2_type_real) then

        call adios2_attribute_data(tmp_attr32, adios2_attr, ierr)
        attr = real(tmp_attr32, real64)

      else if (adios2_attr%type == adios2_type_dp) then

        call adios2_attribute_data(attr, adios2_attr, ierr)

      else
        call error_abort("IO Error: unsuported real type")
      end if

    class default
      call error_abort("Unknown IO process handler type")

    end select

  end subroutine

  !v Read a 1D 32-bit integer array from file
  !
  !  @todo Check if the "mode" can be read from the configuration file
  module subroutine read_array_int32_1D(io_proc, var_name, global_start, count, var, step)
    class(io_process), intent(in) :: io_proc                 !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: var_name                 !< Name of integer array to read
    integer(int64), dimension(1), intent(in) :: global_start !< What global index to start reading from
    integer(int64), dimension(1), intent(in) :: count        !< How many array element to read
    integer(int32), dimension(:), intent(inout) :: var       !< The 1D integer array
    integer(int64), optional, intent(in) :: step             !< The step to read

    type(adios2_variable) :: adios2_var
    integer(int64), dimension(:), allocatable :: tmp_var64
    integer(ccs_int) :: ierr
    integer(int64), parameter :: step_count = 1

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
      call adios2_set_selection(adios2_var, 1, global_start, count, ierr)

      if (present(step)) then
        call adios2_set_step_selection(adios2_var, step, step_count, ierr)
      endif

      if (adios2_var%type == adios2_type_integer8) then

        allocate (tmp_var64(size(var)))

        call downcast_warning()
        call adios2_get(io_proc%engine, adios2_var, tmp_var64, adios2_mode_sync, ierr)
        var = int(tmp_var64, int32)

      else if (adios2_var%type == adios2_type_integer4) then

        call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

      else
        call error_abort("IO Error: unsuported integer type")
      end if

    class default
      call error_abort("Unknown IO process handler type")

    end select

  end subroutine

  !v Read a 1D 64-bit integer array from file
  !
  !  @todo Check if the "mode" can be read from the configuration file
  module subroutine read_array_int64_1D(io_proc, var_name, global_start, count, var, step)
    class(io_process), intent(in) :: io_proc                 !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: var_name                 !< Name of integer array to read
    integer(int64), dimension(1), intent(in) :: global_start !< What global index to start reading from
    integer(int64), dimension(1), intent(in) :: count        !< How many array element to read
    integer(int64), dimension(:), intent(inout) :: var       !< The 1D integer array
    integer(int64), optional, intent(in) :: step             !< The step to read

    type(adios2_variable) :: adios2_var
    integer(int32), dimension(:), allocatable :: tmp_var32
    integer(ccs_int) :: ierr
    integer(int64), parameter :: step_count = 1

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
      call adios2_set_selection(adios2_var, 1, global_start, count, ierr)

      if (present(step)) then
        call adios2_set_step_selection(adios2_var, step, step_count, ierr)
      endif

      if (adios2_var%type == adios2_type_integer4) then

        allocate (tmp_var32(size(var)))

        call adios2_get(io_proc%engine, adios2_var, tmp_var32, adios2_mode_sync, ierr)
        var = int(tmp_var32, int64)

      else if (adios2_var%type == adios2_type_integer8) then

        call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

      else
        call error_abort("IO Error: unsuported integer type")
      end if

    class default
      call error_abort("Unknown IO process handler type")

    end select

  end subroutine

  !v Read a 2D 32-bit integer array from file
  !
  !  @todo Check if the "mode" can be read from the configuration file
  module subroutine read_array_int32_2D(io_proc, var_name, global_start, count, var, step)
    class(io_process), intent(in) :: io_proc                 !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: var_name                 !< Name of integer array to read
    integer(int64), dimension(2), intent(in) :: global_start !< What global index to start reading from
    integer(int64), dimension(2), intent(in) :: count        !< How many array elements to read
    integer(int32), dimension(:, :), intent(inout) :: var    !< The 2D integer array
    integer(int64), optional, intent(in) :: step             !< The step to read

    type(adios2_variable) :: adios2_var
    integer(int64), dimension(:, :), allocatable :: tmp_var64
    integer(ccs_int) :: ierr
    integer(int64), parameter :: step_count = 1

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
      call adios2_set_selection(adios2_var, 2, global_start, count, ierr)

      if (present(step)) then
        call adios2_set_step_selection(adios2_var, step, step_count, ierr)
      endif

      if (adios2_var%type == adios2_type_integer8) then

        allocate (tmp_var64(size(var, dim=1), size(var, dim=2)))

        call downcast_warning()
        call adios2_get(io_proc%engine, adios2_var, tmp_var64, adios2_mode_sync, ierr)
        var = int(tmp_var64, int32)

      else if (adios2_var%type == adios2_type_integer4) then

        call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

      else
        call error_abort("IO Error: unsuported integer type")
      end if

    class default
      call error_abort("Unknown IO process handler type")

    end select

  end subroutine

  !v Read a 2D 64-bit integer array from file
  !
  !  @todo Check if the "mode" can be read from the configuration file
  module subroutine read_array_int64_2D(io_proc, var_name, global_start, count, var, step)
    class(io_process), intent(in) :: io_proc                 !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: var_name                 !< Name of integer array to read
    integer(int64), dimension(2), intent(in) :: global_start !< What global index to start reading from
    integer(int64), dimension(2), intent(in) :: count        !< How many array element to read
    integer(int64), dimension(:, :), intent(inout) :: var    !< The 2D integer array
    integer(int64), optional, intent(in) :: step             !< The step to read

    type(adios2_variable) :: adios2_var
    integer(int32), dimension(:, :), allocatable :: tmp_var32
    integer(ccs_int) :: ierr
    integer(int64), parameter :: step_count = 1

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
      call adios2_set_selection(adios2_var, 2, global_start, count, ierr)

      if (present(step)) then
        call adios2_set_step_selection(adios2_var, step, step_count, ierr)
      endif

      if (adios2_var%type == adios2_type_integer4) then

        allocate (tmp_var32(size(var, dim=1), size(var, dim=2)))

        call adios2_get(io_proc%engine, adios2_var, tmp_var32, adios2_mode_sync, ierr)
        var = int(tmp_var32, int64)

      else if (adios2_var%type == adios2_type_integer8) then

        call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

      else
        call error_abort("IO Error: unsuported integer type")
      end if

    class default
      call error_abort("Unknown IO process handler type")

    end select

  end subroutine

  !v Read a 1D 32-bit real array from file
  !
  !  @todo Check if the "mode" can be read from the configuration file
  module subroutine read_array_real32_1D(io_proc, var_name, global_start, count, var, step)
    class(io_process), intent(in) :: io_proc                 !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: var_name                 !< Name of real array to read
    integer(int64), dimension(1), intent(in) :: global_start !< What global index to start reading from
    integer(int64), dimension(1), intent(in) :: count        !< How many array element to read
    real(real32), dimension(:), intent(inout) :: var         !< The 1D real array
    integer(int64), optional, intent(in) :: step             !< The step to read

    type(adios2_variable) :: adios2_var
    real(real64), dimension(:), allocatable :: tmp_var64
    integer(ccs_int) :: ierr
    integer(int64), parameter :: step_count = 1

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
      call adios2_set_selection(adios2_var, 1, global_start, count, ierr)

      if (present(step)) then
        call adios2_set_step_selection(adios2_var, step, step_count, ierr)
      endif

      if (adios2_var%type == adios2_type_dp) then

        allocate (tmp_var64(size(var)))

        call downcast_warning()
        call adios2_get(io_proc%engine, adios2_var, tmp_var64, adios2_mode_sync, ierr)
        var = real(tmp_var64, real32)

      else if (adios2_var%type == adios2_type_real) then

        call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

      else
        call error_abort("IO Error: unsuported real type")
      end if

    class default
      call error_abort("Unknown IO process handler type")

    end select

  end subroutine

  !v Read a 1D 64-bit real array from file
  !
  !  @todo Check if the "mode" can be read from the configuration file
  module subroutine read_array_real64_1D(io_proc, var_name, global_start, count, var, step)
    class(io_process), intent(in) :: io_proc                 !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: var_name                 !< Name of real array to read
    integer(int64), dimension(1), intent(in) :: global_start !< What global index to start reading from
    integer(int64), dimension(1), intent(in) :: count        !< How many array element to read
    real(real64), dimension(:), intent(inout) :: var         !< The 1D real array
    integer(int64), optional, intent(in) :: step             !< The step to read

    type(adios2_variable) :: adios2_var
    real(real32), dimension(:), allocatable :: tmp_var32
    integer(ccs_int) :: ierr
    integer(int64), parameter :: step_count = 1

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
      call adios2_set_selection(adios2_var, 1, global_start, count, ierr)

      if (present(step)) then
        call adios2_set_step_selection(adios2_var, step, step_count, ierr)
      endif

      if (adios2_var%type == adios2_type_real) then

        allocate (tmp_var32(size(var)))

        call adios2_get(io_proc%engine, adios2_var, tmp_var32, adios2_mode_sync, ierr)
        var = real(tmp_var32, real64)

      else if (adios2_var%type == adios2_type_dp) then

        call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

      else
        call error_abort("IO Error: unsuported real type")
      end if

    class default
      call error_abort("Unknown IO process handler type")

    end select

  end subroutine

  !v Read a 2D 32-bit real array from file
  !
  !  @todo Check if the "mode" can be read from the configuration file
  module subroutine read_array_real32_2D(io_proc, var_name, global_start, count, var, step)
    class(io_process), intent(in) :: io_proc                 !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: var_name                 !< Name of real array to read
    integer(int64), dimension(2), intent(in) :: global_start !< What global index to start reading from
    integer(int64), dimension(2), intent(in) :: count        !< How many array element to read
    real(real32), dimension(:, :), intent(inout) :: var      !< The 2D real array
    integer(int64), optional, intent(in) :: step             !< The step to read

    type(adios2_variable) :: adios2_var
    real(real64), dimension(:, :), allocatable :: tmp_var64
    integer(ccs_int) :: ierr
    integer(int64), parameter :: step_count = 1

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
      call adios2_set_selection(adios2_var, 2, global_start, count, ierr)

      if (present(step)) then
        call adios2_set_step_selection(adios2_var, step, step_count, ierr)
      endif

      if (adios2_var%type == adios2_type_dp) then

        allocate (tmp_var64(size(var, dim=1), size(var, dim=2)))

        call downcast_warning()
        call adios2_get(io_proc%engine, adios2_var, tmp_var64, adios2_mode_sync, ierr)
        var = real(tmp_var64, real32)

      else if (adios2_var%type == adios2_type_real) then

        call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

      else
        call error_abort("IO Error: unsuported real type")
      end if

    class default
      call error_abort("Unknown IO process handler type")

    end select

  end subroutine

  !v Read a 2D 64-bit real array from file
  !
  !  @todo Check if the "mode" can be read from the configuration file
  module subroutine read_array_real64_2D(io_proc, var_name, global_start, count, var, step)
    class(io_process), intent(in) :: io_proc                 !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: var_name                 !< Name of real array to read
    integer(int64), dimension(2), intent(in) :: global_start !< What global index to start reading from
    integer(int64), dimension(2), intent(in) :: count        !< How many array element to read
    real(real64), dimension(:, :), intent(inout) :: var      !< The 2D real array
    integer(int64), optional, intent(in) :: step             !< The step to read

    type(adios2_variable) :: adios2_var
    real(real32), dimension(:, :), allocatable :: tmp_var32
    integer(ccs_int) :: ierr
    integer(int64), parameter :: step_count = 1

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)

      call adios2_set_selection(adios2_var, 2, global_start, count, ierr)

      if (present(step)) then
        call adios2_set_step_selection(adios2_var, step, step_count, ierr)
      endif

      if (adios2_var%type == adios2_type_real) then

        allocate (tmp_var32(size(var, dim=1), size(var, dim=2)))

        call adios2_get(io_proc%engine, adios2_var, tmp_var32, adios2_mode_sync, ierr)
        var = real(tmp_var32, real64)

      else if (adios2_var%type == adios2_type_dp) then

        call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

      else
        call error_abort("IO Error: unsuported real type")
      end if

    class default
      call error_abort("Unknown IO process handler type")

    end select

  end subroutine

  !>  Write a scalar 32-bit integer to file
  module subroutine write_scalar_int32(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc     !< ADIOS2 IO process used for writing
    character(len=*), intent(in) :: attr_name   !< Name of scalar integer to write
    integer(int32), intent(in) :: attr           !< Value of scalar integer

    type(adios2_attribute) :: adios2_attr

    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_define_attribute(adios2_attr, io_proc%io_task, attr_name, attr, ierr)

    class default
      print *, "Unknown IO process handler type"
      stop 1

    end select

  end subroutine

  !>  Write a scalar 64-bit integer from file
  module subroutine write_scalar_int64(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc     !< ADIOS2 IO process used for writing
    character(len=*), intent(in) :: attr_name   !< Name of scalar longinteger to write
    integer(int64), intent(in) :: attr           !< Value of scalar long integer

    type(adios2_attribute) :: adios2_attr

    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_define_attribute(adios2_attr, io_proc%io_task, attr_name, attr, ierr)

    class default
      print *, "Unknown IO process handler type"
      stop 1

    end select

  end subroutine

  !>  Write a scalar 32-bit real from file
  module subroutine write_scalar_real32(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc      !< ADIOS2 IO process used for writing
    character(len=*), intent(in) :: attr_name    !< Name of scalar real to write
    real(real32), intent(in) :: attr             !< Value of scalar real

    type(adios2_attribute) :: adios2_attr

    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_define_attribute(adios2_attr, io_proc%io_task, attr_name, attr, ierr)

    class default
      print *, "Unknown IO process handler type"
      stop 1

    end select

  end subroutine

  !>  Write a scalar 64-bit real from file
  module subroutine write_scalar_real64(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc     !< ADIOS2 IO process used for writing
    character(len=*), intent(in) :: attr_name   !< Name of scalar double precision real to write
    real(real64), intent(in) :: attr             !< Value of scalar double precision real

    type(adios2_attribute) :: adios2_attr

    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_define_attribute(adios2_attr, io_proc%io_task, attr_name, attr, ierr)

    class default
      print *, "Unknown IO process handler type"
      stop 1

    end select

  end subroutine

  !v Write a 1D 32-bit integer array to file
  !
  !  @todo Check if the "mode" can be read from the configuration file
  module subroutine write_array_int32_1D(io_proc, var_name, global_shape, global_start, count, var)
    class(io_process), intent(in) :: io_proc                  !< ADIOS2 IO process used for writing
    character(len=*), intent(in) :: var_name                 !< Name of integer array to write
    integer(int64), dimension(1), intent(in) :: global_shape  !< Global shape of array
    integer(int64), dimension(1), intent(in) :: global_start  !< What global index to start writing from
    integer(int64), dimension(1), intent(in) :: count         !< How many array element to write
    integer(int32), dimension(:), intent(in) :: var           !< The 1D integer array

    type(adios2_variable) :: adios2_var
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_define_variable(adios2_var, io_proc%io_task, var_name, adios2_type_integer4, &
                                  1, global_shape, global_start, count, &
                                  adios2_constant_dims, ierr)
      call adios2_put(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

    class default
      print *, "Unknown IO process handler type"
      stop 1

    end select

  end subroutine

  !v Write a 1D 64-bit integer array from file
  !
  !  @todo Check if the "mode" can be read from the configuration file
  module subroutine write_array_int64_1D(io_proc, var_name, global_shape, global_start, count, var)
    class(io_process), intent(in) :: io_proc                  !< ADIOS2 IO process used for writing
    character(len=*), intent(in) :: var_name                 !< Name of integer array to write
    integer(int64), dimension(1), intent(in) :: global_shape  !< Global shape of array
    integer(int64), dimension(1), intent(in) :: global_start  !< What global index to start writing from
    integer(int64), dimension(1), intent(in) :: count         !< How many array element to write
    integer(int64), dimension(:), intent(in) :: var           !< The 1D integer array

    type(adios2_variable) :: adios2_var
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_define_variable(adios2_var, io_proc%io_task, var_name, adios2_type_integer8, &
                                  1, global_shape, global_start, count, &
                                  adios2_constant_dims, ierr)
      call adios2_put(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

    class default
      print *, "Unknown IO process handler type"
      stop 1

    end select

  end subroutine

  !v Write a 2D 32-bit integer array from file
  !
  !  @todo Check if the "mode" can be read from the configuration file
  module subroutine write_array_int32_2D(io_proc, var_name, global_shape, global_start, count, var)
    class(io_process), intent(in) :: io_proc                  !< ADIOS2 IO process used for writing
    character(len=*), intent(in) :: var_name                 !< Name of integer array to write
    integer(int64), dimension(2), intent(in) :: global_shape  !< Global shape of array
    integer(int64), dimension(2), intent(in) :: global_start  !< What global index to start writing from
    integer(int64), dimension(2), intent(in) :: count         !< How many array elements to write
    integer(int32), dimension(:, :), intent(in) :: var         !< The 2D integer array

    type(adios2_variable) :: adios2_var
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_define_variable(adios2_var, io_proc%io_task, var_name, adios2_type_integer4, &
                                  2, global_shape, global_start, count, &
                                  adios2_constant_dims, ierr)
      call adios2_put(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

    class default
      print *, "Unknown IO process handler type"
      stop 1

    end select

  end subroutine

  !v Write a 2D 64-bit integer array from file
  !
  !  @todo Check if the "mode" can be read from the configuration file
  module subroutine write_array_int64_2D(io_proc, var_name, global_shape, global_start, count, var)
    class(io_process), intent(in) :: io_proc                  !< ADIOS2 IO process used for writing
    character(len=*), intent(in) :: var_name                 !< Name of integer array to write
    integer(int64), dimension(2), intent(in) :: global_shape  !< Global shape of array
    integer(int64), dimension(2), intent(in) :: global_start  !< What global index to start writing from
    integer(int64), dimension(2), intent(in) :: count         !< How many array element to write
    integer(int64), dimension(:, :), intent(in) :: var         !< The 2D integer array

    type(adios2_variable) :: adios2_var
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_define_variable(adios2_var, io_proc%io_task, var_name, adios2_type_integer8, &
                                  2, global_shape, global_start, count, &
                                  adios2_constant_dims, ierr)
      call adios2_put(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

    class default
      print *, "Unknown IO process handler type"
      stop 1

    end select

  end subroutine

  !v Write a 1D 32-bit real array from file
  !
  !  @todo Check if the "mode" can be read from the configuration file
  module subroutine write_array_real32_1D(io_proc, var_name, global_shape, global_start, count, var)
    class(io_process), intent(in) :: io_proc                  !< ADIOS2 IO process used for writing
    character(len=*), intent(in) :: var_name                 !< Name of real array to write
    integer(int64), dimension(1), intent(in) :: global_shape  !< Global shape of array
    integer(int64), dimension(1), intent(in) :: global_start  !< What global index to start writing from
    integer(int64), dimension(1), intent(in) :: count         !< How many array element to write
    real(real32), dimension(:), intent(in) :: var             !< The 1D real array

    type(adios2_variable) :: adios2_var
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_define_variable(adios2_var, io_proc%io_task, var_name, adios2_type_real, &
                                  1, global_shape, global_start, count, &
                                  adios2_constant_dims, ierr)
      call adios2_put(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

    class default
      print *, "Unknown IO process handler type"
      stop 1

    end select

  end subroutine

  !v Write a 1D 64-bit real array from file
  !
  !  @todo Check if the "mode" can be read from the configuration file
  module subroutine write_array_real64_1D(io_proc, var_name, global_shape, global_start, count, var)
    class(io_process), intent(in) :: io_proc                  !< ADIOS2 IO process used for writing
    character(len=*), intent(in) :: var_name                 !< Name of real array to write
    integer(int64), dimension(1), intent(in) :: global_shape  !< Global shape of array
    integer(int64), dimension(1), intent(in) :: global_start  !< What global index to start writing from
    integer(int64), dimension(1), intent(in) :: count         !< How many array element to write
    real(real64), dimension(:), intent(inout) :: var             !< The 1D real array

    type(adios2_variable) :: adios2_var
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)

      if (.not. adios2_var%valid) then
        call adios2_define_variable(adios2_var, io_proc%io_task, var_name, adios2_type_dp, &
                                    1, global_shape, global_start, count, &
                                    adios2_constant_dims, ierr)
      end if
      call adios2_put(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

    class default
      print *, "Unknown IO process handler type"
      stop 1

    end select

  end subroutine

  !v Write a 2D 32-bit real array from file
  !
  !  @todo Check if the "mode" can be read from the configuration file
  module subroutine write_array_real32_2D(io_proc, var_name, global_shape, global_start, count, var)
    class(io_process), intent(in) :: io_proc                  !< ADIOS2 IO process used for writing
    character(len=*), intent(in) :: var_name                 !< Name of real array to write
    integer(int64), dimension(2), intent(in) :: global_shape  !< Global shape of array
    integer(int64), dimension(2), intent(in) :: global_start  !< What global index to start writing from
    integer(int64), dimension(2), intent(in) :: count         !< How many array element to write
    real(real32), dimension(:, :), intent(in) :: var           !< The 2D real array

    type(adios2_variable) :: adios2_var
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_define_variable(adios2_var, io_proc%io_task, var_name, adios2_type_real, &
                                  2, global_shape, global_start, count, &
                                  adios2_constant_dims, ierr)
      call adios2_put(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

    class default
      print *, "Unknown IO process handler type"
      stop 1

    end select

  end subroutine

  !v Write a 2D 64-bit real array from file
  !
  !  @todo Check if the "mode" can be read from the configuration file
  module subroutine write_array_real64_2D(io_proc, var_name, global_shape, global_start, count, var)
    class(io_process), intent(in) :: io_proc                  !< ADIOS2 IO process used for writing
    character(len=*), intent(in) :: var_name                 !< Name of real array to write
    integer(int64), dimension(2), intent(in) :: global_shape  !< Global shape of array
    integer(int64), dimension(2), intent(in) :: global_start  !< What global index to start writing from
    integer(int64), dimension(2), intent(in) :: count         !< How many array element to write
    real(real64), dimension(:, :), intent(in) :: var           !< The 2D real array

    type(adios2_variable) :: adios2_var
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_define_variable(adios2_var, io_proc%io_task, var_name, adios2_type_dp, &
                                  2, global_shape, global_start, count, &
                                  adios2_constant_dims, ierr)
      call adios2_put(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

    class default
      print *, "Unknown IO process handler type"
      stop 1

    end select

  end subroutine

  !>  Print out downcast warning
  subroutine downcast_warning()
    print *, "===> IO Warning:"
    print *, "===> Downcasting from 64-bit to 32-bit, possible loss of precision."
  end subroutine

end submodule
