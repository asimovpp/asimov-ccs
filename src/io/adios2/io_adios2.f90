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
  use kinds, only: ccs_int, ccs_real, ccs_long

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
  module subroutine read_array_int32_1D(io_proc, var_name, global_start, count, var)
    class(io_process), intent(in) :: io_proc                 !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: var_name                 !< Name of integer array to read
    integer(int64), dimension(1), intent(in) :: global_start !< What global index to start reading from
    integer(int64), dimension(1), intent(in) :: count        !< How many array element to read
    integer(int32), dimension(:), intent(inout) :: var       !< The 1D integer array

    type(adios2_variable) :: adios2_var
    integer(int64), dimension(:), allocatable :: tmp_var64
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
      call adios2_set_selection(adios2_var, 1, global_start, count, ierr)

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
  module subroutine read_array_int64_1D(io_proc, var_name, global_start, count, var)
    class(io_process), intent(in) :: io_proc                 !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: var_name                 !< Name of integer array to read
    integer(int64), dimension(1), intent(in) :: global_start !< What global index to start reading from
    integer(int64), dimension(1), intent(in) :: count        !< How many array element to read
    integer(int64), dimension(:), intent(inout) :: var       !< The 1D integer array

    type(adios2_variable) :: adios2_var
    integer(int32), dimension(:), allocatable :: tmp_var32
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
      call adios2_set_selection(adios2_var, 1, global_start, count, ierr)

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
  module subroutine read_array_int32_2D(io_proc, var_name, global_start, count, var)
    class(io_process), intent(in) :: io_proc                 !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: var_name                 !< Name of integer array to read
    integer(int64), dimension(2), intent(in) :: global_start !< What global index to start reading from
    integer(int64), dimension(2), intent(in) :: count        !< How many array elements to read
    integer(int32), dimension(:, :), intent(inout) :: var    !< The 2D integer array

    type(adios2_variable) :: adios2_var
    integer(int64), dimension(:, :), allocatable :: tmp_var64
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
      call adios2_set_selection(adios2_var, 2, global_start, count, ierr)

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
  module subroutine read_array_int64_2D(io_proc, var_name, global_start, count, var)
    class(io_process), intent(in) :: io_proc                 !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: var_name                 !< Name of integer array to read
    integer(int64), dimension(2), intent(in) :: global_start !< What global index to start reading from
    integer(int64), dimension(2), intent(in) :: count        !< How many array element to read
    integer(int64), dimension(:, :), intent(inout) :: var    !< The 2D integer array

    type(adios2_variable) :: adios2_var
    integer(int32), dimension(:, :), allocatable :: tmp_var32
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
      call adios2_set_selection(adios2_var, 2, global_start, count, ierr)

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
  module subroutine read_array_real32_1D(io_proc, var_name, global_start, count, var)
    class(io_process), intent(in) :: io_proc                 !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: var_name                 !< Name of real array to read
    integer(int64), dimension(1), intent(in) :: global_start !< What global index to start reading from
    integer(int64), dimension(1), intent(in) :: count        !< How many array element to read
    real(real32), dimension(:), intent(inout) :: var         !< The 1D real array

    type(adios2_variable) :: adios2_var
    real(real64), dimension(:), allocatable :: tmp_var64
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
      call adios2_set_selection(adios2_var, 1, global_start, count, ierr)

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
  module subroutine read_array_real64_1D(io_proc, var_name, global_start, count, var)
    class(io_process), intent(in) :: io_proc                 !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: var_name                 !< Name of real array to read
    integer(int64), dimension(1), intent(in) :: global_start !< What global index to start reading from
    integer(int64), dimension(1), intent(in) :: count        !< How many array element to read
    real(real64), dimension(:), intent(inout) :: var         !< The 1D real array

    type(adios2_variable) :: adios2_var
    real(real32), dimension(:), allocatable :: tmp_var32
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
      call adios2_set_selection(adios2_var, 1, global_start, count, ierr)

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
  module subroutine read_array_real32_2D(io_proc, var_name, global_start, count, var)
    class(io_process), intent(in) :: io_proc                 !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: var_name                 !< Name of real array to read
    integer(int64), dimension(2), intent(in) :: global_start !< What global index to start reading from
    integer(int64), dimension(2), intent(in) :: count        !< How many array element to read
    real(real32), dimension(:, :), intent(inout) :: var      !< The 2D real array

    type(adios2_variable) :: adios2_var
    real(real64), dimension(:, :), allocatable :: tmp_var64
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
      call adios2_set_selection(adios2_var, 2, global_start, count, ierr)

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
  module subroutine read_array_real64_2D(io_proc, var_name, global_start, count, var)
    class(io_process), intent(in) :: io_proc                 !< ADIOS2 IO process used for reading
    character(len=*), intent(in) :: var_name                 !< Name of real array to read
    integer(int64), dimension(2), intent(in) :: global_start !< What global index to start reading from
    integer(int64), dimension(2), intent(in) :: count        !< How many array element to read
    real(real64), dimension(:, :), intent(inout) :: var      !< The 2D real array

    type(adios2_variable) :: adios2_var
    real(real32), dimension(:, :), allocatable :: tmp_var32
    integer(ccs_int) :: ierr

    select type (io_proc)
    type is (adios2_io_process)

      call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)

      call adios2_set_selection(adios2_var, 2, global_start, count, ierr)
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
    character (len=*), intent(in) :: attr_name   !< Name of scalar integer to write
    integer(int32), intent(in) :: attr           !< Value of scalar integer

    type(adios2_attribute) :: adios2_attr

    integer(ccs_int) :: ierr

    select type(io_proc)
      type is(adios2_io_process)

        call adios2_define_attribute(adios2_attr, io_proc%io_task, attr_name, attr, ierr)

      class default
        print*,"Unknown IO process handler type"
        stop 1

      end select

  end subroutine

  !>  Write a scalar 64-bit integer from file
  module subroutine write_scalar_int64(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc     !< ADIOS2 IO process used for writing
    character (len=*), intent(in) :: attr_name   !< Name of scalar longinteger to write
    integer(int64), intent(in) :: attr           !< Value of scalar long integer

    type(adios2_attribute) :: adios2_attr

    integer(ccs_int) :: ierr

    select type(io_proc)
      type is(adios2_io_process)

        call adios2_define_attribute(adios2_attr, io_proc%io_task, attr_name, attr, ierr)

      class default
        print*,"Unknown IO process handler type"
        stop 1

      end select

  end subroutine

  !>  Write a scalar 32-bit real from file
  module subroutine write_scalar_real32(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc      !< ADIOS2 IO process used for writing
    character (len=*), intent(in) :: attr_name    !< Name of scalar real to write
    real(real32), intent(in) :: attr             !< Value of scalar real

    type(adios2_attribute) :: adios2_attr

    integer(ccs_int) :: ierr

    select type(io_proc)
      type is(adios2_io_process)

        call adios2_define_attribute(adios2_attr, io_proc%io_task, attr_name, attr, ierr)

      class default
        print*,"Unknown IO process handler type"
        stop 1

      end select

  end subroutine

  !>  Write a scalar 64-bit real from file
  module subroutine write_scalar_real64(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc     !< ADIOS2 IO process used for writing
    character (len=*), intent(in) :: attr_name   !< Name of scalar double precision real to write
    real(real64), intent(in) :: attr             !< Value of scalar double precision real

    type(adios2_attribute) :: adios2_attr

    integer(ccs_int) :: ierr

    select type(io_proc)
      type is(adios2_io_process)

        call adios2_define_attribute(adios2_attr, io_proc%io_task, attr_name, attr, ierr)
  
      class default
        print*,"Unknown IO process handler type"
        stop 1

      end select

    end subroutine
  

    !>  Write a 1D 32-bit integer array to file
    !
    !> @todo Check if the "mode" can be read from the configuration file
    module subroutine write_array_int32_1D(io_proc, var_name, global_shape, global_start, count, var)
      class(io_process), intent(in) :: io_proc                  !< ADIOS2 IO process used for writing
      character (len=*), intent(in) :: var_name                 !< Name of integer array to write
      integer(int64), dimension(1), intent(in) :: global_shape  !< Global shape of array
      integer(int64), dimension(1), intent(in) :: global_start  !< What global index to start writing from
      integer(int64), dimension(1), intent(in) :: count         !< How many array element to write
      integer(int32), dimension(:), intent(in) :: var           !< The 1D integer array

      type(adios2_variable):: adios2_var
      integer(ccs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_define_variable(adios2_var, io_proc%io_task, var_name, adios2_type_integer4, &
                                      1, global_shape, global_start, count, &
                                      adios2_constant_dims, ierr)
          call adios2_put(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)
  
      class default
        print*,"Unknown IO process handler type"
        stop 1

      end select

    end subroutine

    !>  Write a 1D 64-bit integer array from file
    !
    !> @todo Check if the "mode" can be read from the configuration file
    module subroutine write_array_int64_1D(io_proc, var_name, global_shape, global_start, count, var)
      class(io_process), intent(in) :: io_proc                  !< ADIOS2 IO process used for writing
      character (len=*), intent(in) :: var_name                 !< Name of integer array to write
      integer(int64), dimension(1), intent(in) :: global_shape  !< Global shape of array
      integer(int64), dimension(1), intent(in) :: global_start  !< What global index to start writing from
      integer(int64), dimension(1), intent(in) :: count         !< How many array element to write
      integer(int64), dimension(:), intent(in) :: var           !< The 1D integer array    

      type(adios2_variable):: adios2_var
      integer(ccs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_define_variable(adios2_var, io_proc%io_task, var_name, adios2_type_integer8, &
                                      1, global_shape, global_start, count, &
                                      adios2_constant_dims, ierr)
          call adios2_put(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)
  
      class default
        print*,"Unknown IO process handler type"
        stop 1

      end select

    end subroutine

    !>  Write a 2D 32-bit integer array from file
    !
    !> @todo Check if the "mode" can be read from the configuration file
    module subroutine write_array_int32_2D(io_proc, var_name, global_shape, global_start, count, var)
      class(io_process), intent(in) :: io_proc                  !< ADIOS2 IO process used for writing
      character (len=*), intent(in) :: var_name                 !< Name of integer array to write
      integer(int64), dimension(2), intent(in) :: global_shape  !< Global shape of array
      integer(int64), dimension(2), intent(in) :: global_start  !< What global index to start writing from
      integer(int64), dimension(2), intent(in) :: count         !< How many array elements to write
      integer(int32), dimension(:,:), intent(in) :: var         !< The 2D integer array

      type(adios2_variable):: adios2_var
      integer(ccs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_define_variable(adios2_var, io_proc%io_task, var_name, adios2_type_integer4, &
                                      2, global_shape, global_start, count, &
                                      adios2_constant_dims, ierr)
          call adios2_put(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

      class default
        print*,"Unknown IO process handler type"
        stop 1

      end select

    end subroutine

    !>  Write a 2D 64-bit integer array from file
    !
    !> @todo Check if the "mode" can be read from the configuration file
    module subroutine write_array_int64_2D(io_proc, var_name, global_shape, global_start, count, var)
      class(io_process), intent(in) :: io_proc                  !< ADIOS2 IO process used for writing
      character (len=*), intent(in) :: var_name                 !< Name of integer array to write
      integer(int64), dimension(2), intent(in) :: global_shape  !< Global shape of array
      integer(int64), dimension(2), intent(in) :: global_start  !< What global index to start writing from
      integer(int64), dimension(2), intent(in) :: count         !< How many array element to write
      integer(int64), dimension(:,:), intent(in) :: var         !< The 2D integer array

      type(adios2_variable):: adios2_var
      integer(ccs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_define_variable(adios2_var, io_proc%io_task, var_name, adios2_type_integer8, &
                                      2, global_shape, global_start, count, &
                                      adios2_constant_dims, ierr)
          call adios2_put(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

      class default
        print*,"Unknown IO process handler type"
        stop 1

      end select

    end subroutine

    !>  Write a 1D 32-bit real array from file
    !
    !> @todo Check if the "mode" can be read from the configuration file
    module subroutine write_array_real32_1D(io_proc, var_name, global_shape, global_start, count, var)
      class(io_process), intent(in) :: io_proc                  !< ADIOS2 IO process used for writing
      character (len=*), intent(in) :: var_name                 !< Name of real array to write
      integer(int64), dimension(1), intent(in) :: global_shape  !< Global shape of array
      integer(int64), dimension(1), intent(in) :: global_start  !< What global index to start writing from
      integer(int64), dimension(1), intent(in) :: count         !< How many array element to write
      real(real32), dimension(:), intent(in) :: var             !< The 1D real array

      type(adios2_variable):: adios2_var
      integer(ccs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_define_variable(adios2_var, io_proc%io_task, var_name, adios2_type_real, &
                                      1, global_shape, global_start, count, &
                                      adios2_constant_dims, ierr)
          call adios2_put(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

        class default
          print*,"Unknown IO process handler type"
          stop 1

      end select

    end subroutine

    !>  Write a 1D 64-bit real array from file
    !
    !> @todo Check if the "mode" can be read from the configuration file
    module subroutine write_array_real64_1D(io_proc, var_name, global_shape, global_start, count, var)
      class(io_process), intent(in) :: io_proc                  !< ADIOS2 IO process used for writing
      character (len=*), intent(in) :: var_name                 !< Name of real array to write
      integer(int64), dimension(1), intent(in) :: global_shape  !< Global shape of array
      integer(int64), dimension(1), intent(in) :: global_start  !< What global index to start writing from
      integer(int64), dimension(1), intent(in) :: count         !< How many array element to write
      real(real64), dimension(:), intent(inout) :: var             !< The 1D real array

      type(adios2_variable):: adios2_var
      integer(ccs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)

          if (.not. adios2_var%valid) then
            call adios2_define_variable(adios2_var, io_proc%io_task, var_name, adios2_type_dp, &
                                      1, global_shape, global_start, count, &
                                      adios2_constant_dims, ierr)
          endif
          call adios2_put(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)
          
        class default
          print*,"Unknown IO process handler type"
          stop 1

      end select

    end subroutine

    !>  Write a 2D 32-bit real array from file
    !
    !> @todo Check if the "mode" can be read from the configuration file
    module subroutine write_array_real32_2D(io_proc, var_name, global_shape, global_start, count, var)
      class(io_process), intent(in) :: io_proc                  !< ADIOS2 IO process used for writing
      character (len=*), intent(in) :: var_name                 !< Name of real array to write
      integer(int64), dimension(2), intent(in) :: global_shape  !< Global shape of array
      integer(int64), dimension(2), intent(in) :: global_start  !< What global index to start writing from
      integer(int64), dimension(2), intent(in) :: count         !< How many array element to write
      real(real32), dimension(:,:), intent(in) :: var           !< The 2D real array

      type(adios2_variable):: adios2_var
      integer(ccs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_define_variable(adios2_var, io_proc%io_task, var_name, adios2_type_real, &
                                      2, global_shape, global_start, count, &
                                      adios2_constant_dims, ierr)
          call adios2_put(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

      class default
        print*,"Unknown IO process handler type"
        stop 1

      end select

    end subroutine

    !>  Write a 2D 64-bit real array from file
    !
    !> @todo Check if the "mode" can be read from the configuration file
    module subroutine write_array_real64_2D(io_proc, var_name, global_shape, global_start, count, var)
      class(io_process), intent(in) :: io_proc                  !< ADIOS2 IO process used for writing
      character (len=*), intent(in) :: var_name                 !< Name of real array to write
      integer(int64), dimension(2), intent(in) :: global_shape  !< Global shape of array
      integer(int64), dimension(2), intent(in) :: global_start  !< What global index to start writing from
      integer(int64), dimension(2), intent(in) :: count         !< How many array element to write
      real(real64), dimension(:,:), intent(in) :: var           !< The 2D real array

      type(adios2_variable):: adios2_var
      integer(ccs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_define_variable(adios2_var, io_proc%io_task, var_name, adios2_type_dp, &
                                      2, global_shape, global_start, count, &
                                      adios2_constant_dims, ierr)
          call adios2_put(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

      class default
        print*,"Unknown IO process handler type"
        stop 1
        
      end select

    end subroutine

    !>  Print out downcast warning
    subroutine downcast_warning()
      print*,"===> IO Warning:"
      print*,"===> Downcasting from 64-bit to 32-bit, possible loss of precision."
    end subroutine

    subroutine begin_step(io_proc)
      class(io_process), intent(inout) :: io_proc

      integer(ccs_int) :: ierr

      select type (io_proc)
      type is (adios2_io_process)

        call adios2_begin_step(io_proc%engine, ierr)

      class default
        call error_abort("Unknown IO process handler type")

      end select
    end subroutine

    subroutine end_step(io_proc)
      class(io_process), intent(inout) :: io_proc

      integer(ccs_int) :: ierr

      select type (io_proc)
      type is (adios2_io_process)

        call adios2_end_step(io_proc%engine, ierr)

      class default
        call error_abort("Unknown IO process handler type")

      end select
    end subroutine

    module subroutine write_solution(par_env, case_name, step, maxstep, dt, mesh, cps, u, v, w, p)

      use kinds, only: ccs_long
      use constants, only: ndim, adiosconfig
      use vec, only : get_vector_data, restore_vector_data
  
      ! Arguments
      class(parallel_environment), allocatable, target, intent(in) :: par_env
      character(len=:), allocatable, intent(in) :: case_name
      integer(ccs_int), intent(in) :: step
      integer(ccs_int), intent(in) :: maxstep
      real(ccs_real), intent(in) :: dt
      type(ccs_mesh), intent(in) :: mesh
      integer(ccs_int), intent(in) :: cps
      class(field), intent(inout) :: u, v, w, p
  
      ! Local variables
      character(len=:), allocatable :: sol_file
      character(len=:), allocatable :: geo_file
      character(len=:), allocatable :: adios2_file
      character(len=:), allocatable :: xdmf_file
  
      class(io_environment), allocatable, save :: io_env
      class(io_process), allocatable, save :: sol_writer
  
      integer(ccs_long), dimension(1) :: sel_shape
      integer(ccs_long), dimension(1) :: sel_start
      integer(ccs_long), dimension(1) :: sel_count
  
      integer(ccs_long), dimension(2) :: sel2_shape
      integer(ccs_long), dimension(2) :: sel2_start
      integer(ccs_long), dimension(2) :: sel2_count
  
      real(ccs_real), dimension(:), pointer :: data
  
      integer(ccs_int), save :: ioxdmf

      integer(ccs_int) :: ierr
      integer(ccs_int), save :: step_counter = 0
  
      geo_file = case_name//'.geo'
      sol_file = case_name//'.sol.h5'
      adios2_file = case_name//adiosconfig
      xdmf_file = case_name//'.sol.xmf'
  
      if (step == 1) then
        call initialise_io(par_env, adios2_file, io_env)
        call configure_io(io_env, "sol_writer", sol_writer)
        call open_file(sol_file, "write", sol_writer)
      endif
  
      ! 1D data
      sel_shape(1) = mesh%topo%global_num_cells
      sel_start(1) = mesh%topo%global_indices(1) - 1
      sel_count(1) = mesh%topo%local_num_cells
  
      ! 2D data
      sel2_shape(1) = ndim
      sel2_shape(2) = mesh%topo%global_num_cells
      sel2_start(1) = 0
      sel2_start(2) = mesh%topo%global_indices(1) - 1
      sel2_count(1) = ndim
      sel2_count(2) = mesh%topo%local_num_cells
  
      ! Begin step
      call begin_step(sol_writer)
  
      ! Write u-velocity
      call get_vector_data(u%values, data)
      call write_array_real64_1d(sol_writer, "/u", sel_shape, sel_start, sel_count, data)
      call restore_vector_data(u%values, data)
  
      ! Write v-velocity
      call get_vector_data(v%values, data)
      call write_array_real64_1d(sol_writer, "/v", sel_shape, sel_start, sel_count, data)
      call restore_vector_data(v%values, data)

      ! Write w-velocity
      call get_vector_data(w%values, data)
      call write_array_real64_1d(sol_writer, "/w", sel_shape, sel_start, sel_count, data)
      call restore_vector_data(w%values, data)

      ! Write pressure
      call get_vector_data(p%values, data)
      call write_array_real64_1d(sol_writer, "/p", sel_shape, sel_start, sel_count, data)
      call restore_vector_data(p%values, data)

      ! Write gradients
      call get_vector_data(u%x_gradients, data)
      call write_array_real64_1d(sol_writer, "/dudx", sel_shape, sel_start, sel_count, data)
      call restore_vector_data(u%x_gradients, data)

      call get_vector_data(u%y_gradients, data)
      call write_array_real64_1d(sol_writer, "/dudy", sel_shape, sel_start, sel_count, data)
      call restore_vector_data(u%y_gradients, data)

      call get_vector_data(u%z_gradients, data)
      call write_array_real64_1d(sol_writer, "/dudz", sel_shape, sel_start, sel_count, data)
      call restore_vector_data(u%z_gradients, data)

      call get_vector_data(v%x_gradients, data)
      call write_array_real64_1d(sol_writer, "/dvdx", sel_shape, sel_start, sel_count, data)
      call restore_vector_data(v%x_gradients, data)

      call get_vector_data(v%y_gradients, data)
      call write_array_real64_1d(sol_writer, "/dvdy", sel_shape, sel_start, sel_count, data)
      call restore_vector_data(v%y_gradients, data)

      call get_vector_data(v%z_gradients, data)
      call write_array_real64_1d(sol_writer, "/dvdz", sel_shape, sel_start, sel_count, data)
      call restore_vector_data(v%z_gradients, data)

      call get_vector_data(w%x_gradients, data)
      call write_array_real64_1d(sol_writer, "/dwdx", sel_shape, sel_start, sel_count, data)
      call restore_vector_data(w%x_gradients, data)

      call get_vector_data(w%y_gradients, data)
      call write_array_real64_1d(sol_writer, "/dwdy", sel_shape, sel_start, sel_count, data)
      call restore_vector_data(w%y_gradients, data)

      call get_vector_data(w%z_gradients, data)
      call write_array_real64_1d(sol_writer, "/dwdz", sel_shape, sel_start, sel_count, data)
      call restore_vector_data(w%z_gradients, data)

      call get_vector_data(p%x_gradients, data)
      call write_array_real64_1d(sol_writer, "/dpdx", sel_shape, sel_start, sel_count, data)
      call restore_vector_data(p%x_gradients, data)

      call get_vector_data(p%y_gradients, data)
      call write_array_real64_1d(sol_writer, "/dpdy", sel_shape, sel_start, sel_count, data)
      call restore_vector_data(p%y_gradients, data)

      call get_vector_data(p%z_gradients, data)
      call write_array_real64_1d(sol_writer, "/dpdz", sel_shape, sel_start, sel_count, data)
      call restore_vector_data(p%z_gradients, data)

      ! End step
      call end_step(sol_writer)

      if (step == maxstep) then
        ! Close the file and ADIOS2 engine
        call close_file(sol_writer)
  
        ! Finalise the ADIOS2 IO environment
        call cleanup_io(io_env)
      endif
  
      ! Write XML file
      if (step == 1) then
        if (par_env%proc_id == par_env%root) then
          ! Open file
          open(newunit=ioxdmf, file=xdmf_file, status='unknown')
  
          ! Write file contents
          write(ioxdmf, '(a)')        '<?xml version="1.0"?>'
          write(ioxdmf, '(a)')        '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
          write(ioxdmf, '(a)')        '<Xdmf Version="2.0">'
          write(ioxdmf, '(a)')        '  <Domain>'
          write(ioxdmf, '(a)')        '    <Grid Name="Unsteady" GridType="Collection" CollectionType="Temporal">'
        endif
      endif

      associate (ncel => mesh%topo%global_num_cells, &
                 nvrt => mesh%topo%global_num_vertices)

        write(ioxdmf, '(a)')            '      <Grid Name="Mesh">'
        write(ioxdmf, '(a,f10.7,a)')    '        <Time Value = "',step*dt,'" />'

        ! Topology
        if (mesh%topo%vert_per_cell == 4) then
          write(ioxdmf, '(a,i0,a)')   '        <Topology Type="Quadrilateral" NumberOfElements="',ncel,'" BaseOffset="1">'
        else
          write(ioxdmf, '(a,i0,a)')   '        <Topology Type="Hexahedron" NumberOfElements="',ncel,'" BaseOffset="1">'
        endif
        write(ioxdmf, '(a,i0,1x,i0,3(a))') '          <DataItem Dimensions="',ncel,mesh%topo%vert_per_cell,'" Format="HDF">',trim(geo_file),':/Step0/cell/vertices</DataItem>'          
        write(ioxdmf, '(a)')            '        </Topology>'

        ! Geometry
        write(ioxdmf, '(a)')            '        <Geometry Type="XYZ">'
        write(ioxdmf, '(a,i0,1x,i0,3(a))') '          <DataItem Dimensions="',nvrt,ndim,'" Format="HDF">',trim(geo_file),':/Step0/vert</DataItem>'        
        write(ioxdmf, '(a)')            '        </Geometry>'

        ! Velocity vector
        write(ioxdmf, '(a)')            '        <Attribute Name="velocity" AttributeType="Vector" Center="Cell">'
        write(ioxdmf, '(a,i0,1x,i0,a,a)') '          <DataItem Dimensions="',ncel,ndim,'" ItemType="Function"', &
                                                     ' Function="JOIN($0, $1, $2)">'
        write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/u</DataItem>'
        write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/v</DataItem>'
        write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/w</DataItem>'
        write(ioxdmf, '(a)')            '          </DataItem>'
        write(ioxdmf, '(a)')            '        </Attribute>'

        ! Pressure
        write(ioxdmf, '(a)')            '        <Attribute Name="pressure" AttributeType="Scalar" Center="Cell">'
        write(ioxdmf, '(a,i0,3(a),i0,a)') '          <DataItem Dimensions="',ncel,'" Format="HDF">',trim(sol_file),':/Step',step_counter,'/p</DataItem>'          
        write(ioxdmf, '(a)')            '        </Attribute>'

        ! Kinetic Energy
        write(ioxdmf, '(a)')            '        <Attribute Name="kinetic energy" AttributeType="Scalar" Center="Cell">'
        write(ioxdmf, '(a,i0,a,a)')     '          <DataItem Dimensions="',ncel,'" ItemType="Function"', &
                                                   ' Function="0.5 * ($0*$0 + $1*$1 + $2*$2)">'
        write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/u</DataItem>'
        write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/v</DataItem>'
        write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/w</DataItem>'
        write(ioxdmf, '(a)')            '          </DataItem>'
        write(ioxdmf, '(a)')            '        </Attribute>'

        ! Enstrophy
        write(ioxdmf, '(a)')            '        <Attribute Name="enstrophy" AttributeType="Scalar" Center="Cell">'
        write(ioxdmf, '(a,i0,a,a)')     '          <DataItem Dimensions="',ncel,'" ItemType="Function"', &
                                                   ' Function="0.5 * (($5-$3)*($5-$3) + ($1-$4)*($1-$4) + ($2-$0)*($2-$0))">'
        write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/dudy</DataItem>'
        write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/dudz</DataItem>'
        write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/dvdx</DataItem>'
        write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/dvdz</DataItem>'
        write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/dwdx</DataItem>'
        write(ioxdmf, '(a,i0,3(a),i0,a)') '            <DataItem Format="HDF" Dimensions="',ncel,'">',trim(sol_file),':/Step',step_counter,'/dwdy</DataItem>'
        write(ioxdmf, '(a)')            '          </DataItem>'
        write(ioxdmf, '(a)')            '        </Attribute>'

        write(ioxdmf, '(a)')            '      </Grid>'
      
      end associate
  
      ! Close file
      if (step == maxstep) then
        write(ioxdmf, '(a)')          '    </Grid>'
        write(ioxdmf, '(a)')          '  </Domain>'
        write(ioxdmf, '(a)')          '</Xdmf>'
        close(ioxdmf)
      endif

      step_counter = step_counter + 1
  
    end subroutine

end submodule