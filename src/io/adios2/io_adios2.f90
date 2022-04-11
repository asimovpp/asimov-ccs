!> @brief Submodule file io_adios2.smod
!
!> @build mpi adios2
!
!> @details Implementation (using MPI and ADIOS@) of parallel 
!!          IO functionality
submodule (io) io_adios2

  use adios2
  use adios2_types, only: adios2_io_process

  implicit none

  contains

  !> @brief Read a scalar 32-bit integer from file
  !
  !> param[in]  io_proc   : ADIOS2 IO process used for reading
  !> param[in]  attr_name : Name of scalar integer to read
  !> param[out] attr      : Value of scalar integer
  module subroutine read_scalar_int32(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: attr_name
    integer(int32), intent(out) :: attr

    type(adios2_attribute) :: adios2_attr

    integer(int64) :: tmp_attr64
    integer(accs_int) :: ierr

    select type(io_proc)
      type is(adios2_io_process)

        call adios2_inquire_attribute(adios2_attr, io_proc%io_task, attr_name, ierr)

        if (adios2_attr%type == adios2_type_integer8) then

          call downcast_warning()
          call adios2_attribute_data(tmp_attr64, adios2_attr, ierr)
          attr = int(tmp_attr64, int32)

        else if (adios2_attr%type == adios2_type_integer4) then

          call adios2_attribute_data(attr, adios2_attr, ierr)

        else
          print*,"===> IO Error: unsuported integer type"
          stop 1
        end if

      class default
        print*,"Unknown IO process handler type"
        stop 1
      
      end select

  end subroutine

  !> @brief Read a scalar 64-bit integer from file
  !
  !> param[in]  io_proc   : ADIOS2 IO process used for reading
  !> param[in]  attr_name : Name of scalar longinteger to read
  !> param[out] attr      : Value of scalar long integer
  module subroutine read_scalar_int64(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: attr_name
    integer(int64), intent(out) :: attr

    type(adios2_attribute) :: adios2_attr

    integer(int32) :: tmp_attr32
    integer(accs_int) :: ierr

    select type(io_proc)
      type is(adios2_io_process)

        call adios2_inquire_attribute(adios2_attr, io_proc%io_task, attr_name, ierr)

        if (adios2_attr%type == adios2_type_integer4) then

          call adios2_attribute_data(tmp_attr32, adios2_attr, ierr)
          attr = int(tmp_attr32, int64)

        else if (adios2_attr%type == adios2_type_integer8) then

          call adios2_attribute_data(attr, adios2_attr, ierr)

        else
          print*,"===> IO Error: unsuported integer type"
          stop 1
        end if

      class default
        print*,"Unknown IO process handler type"
        stop 1

      end select

  end subroutine

  !> @brief Read a scalar 32-bit real from file
  !
  !> param[in]  io_proc   : ADIOS2 IO process used for reading
  !> param[in]  attr_name : Name of scalar real to read
  !> param[out] attr      : Value of scalar real
  module subroutine read_scalar_real32(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: attr_name
    real(real32), intent(out) :: attr

    type(adios2_attribute) :: adios2_attr

    real(real64) :: tmp_attr64
    integer(accs_int) :: ierr

    select type(io_proc)
      type is(adios2_io_process)

        call adios2_inquire_attribute(adios2_attr, io_proc%io_task, attr_name, ierr)
        if (adios2_attr%type == adios2_type_dp) then
 
          call downcast_warning()
          call adios2_attribute_data(tmp_attr64, adios2_attr, ierr)
          attr = real(tmp_attr64, real32)

        else if (adios2_attr%type == adios2_type_real) then
 
          call adios2_attribute_data(attr, adios2_attr, ierr)
 
        else
          print*,"===> IO Error: unsuported real type"
          stop 1
        end if

      class default
        print*,"Unknown IO process handler type"
        stop 1

      end select

    end subroutine

  !> @brief Read a scalar 64-bit real from file
  !
  !> param[in]  io_proc   : ADIOS2 IO process used for reading
  !> param[in]  attr_name : Name of scalar double precision real to read
  !> param[out] attr      : Value of scalar double precision real
  module subroutine read_scalar_real64(io_proc, attr_name, attr)
    class(io_process), intent(in) :: io_proc
    character (len=*), intent(in) :: attr_name
    real(real64), intent(out) :: attr

    type(adios2_attribute) :: adios2_attr

    real(real32) :: tmp_attr32
    integer(accs_int) :: ierr

    select type(io_proc)
      type is(adios2_io_process)

        call adios2_inquire_attribute(adios2_attr, io_proc%io_task, attr_name, ierr)
        if (adios2_attr%type == adios2_type_real) then

          call adios2_attribute_data(tmp_attr32, adios2_attr, ierr)
          attr = real(tmp_attr32, real64)

        else if (adios2_attr%type == adios2_type_dp) then
  
          call adios2_attribute_data(attr, adios2_attr, ierr)
  
        else
          print*,"===> IO Error: unsuported real type"
          stop 1
        end if
  
      class default
        print*,"Unknown IO process handler type"
        stop 1

      end select

    end subroutine
  

    !> @brief Read a 1D 32-bit integer array from file
    !
    !> @todo Check if the "mode" can be read from the configuration file
    !
    !> param[in]    io_proc  : ADIOS2 IO process used for reading
    !> param[in]    var_name : Name of integer array to read
    !> param[in]    start    : What global index to start reading from
    !> param[in]    count    : How many array element to read
    !> param[input] var      : The 1D integer array
    module subroutine read_array_int32_1D(io_proc, var_name, start, count, var)
      class(io_process), intent(in) :: io_proc
      character (len=*), intent(in) :: var_name
      integer(int64), dimension(1), intent(in) :: start
      integer(int64), dimension(1), intent(in) :: count
      integer(int32), dimension(:), intent(inout) :: var

      type(adios2_variable):: adios2_var
      integer(int64), dimension(:), allocatable :: tmp_var64
      integer(accs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
          call adios2_set_selection(adios2_var, 1, start, count, ierr)

          if (adios2_var%type == adios2_type_integer8) then

            allocate(tmp_var64(size(var)))

            call downcast_warning()
            call adios2_get(io_proc%engine, adios2_var, tmp_var64, adios2_mode_sync, ierr)
            var = int(tmp_var64, int32)
  
          else if (adios2_var%type == adios2_type_integer4) then
  
            call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)
  
          else
            print*,"===> IO Error: unsuported integer type"
            stop 1
          end if
  
      class default
        print*,"Unknown IO process handler type"
        stop 1

      end select

    end subroutine

    !> @brief Read a 1D 64-bit integer array from file
    !
    !> @todo Check if the "mode" can be read from the configuration file
    !
    !> param[in]    io_proc  : ADIOS2 IO process used for reading
    !> param[in]    var_name : Name of integer array to read
    !> param[in]    start    : What global index to start reading from
    !> param[in]    count    : How many array element to read
    !> param[input] var      : The 1D integer array    
    module subroutine read_array_int64_1D(io_proc, var_name, start, count, var)
      class(io_process), intent(in) :: io_proc
      character (len=*), intent(in) :: var_name
      integer(int64), dimension(1), intent(in) :: start
      integer(int64), dimension(1), intent(in) :: count
      integer(int64), dimension(:), intent(inout) :: var

      type(adios2_variable):: adios2_var
      integer(int32), dimension(:), allocatable :: tmp_var32
      integer(accs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
          call adios2_set_selection(adios2_var, 1, start, count, ierr)

          if (adios2_var%type == adios2_type_integer4) then

            allocate(tmp_var32(size(var)))

            call adios2_get(io_proc%engine, adios2_var, tmp_var32, adios2_mode_sync, ierr)
            var = int(tmp_var32, int64)
  
          else if (adios2_var%type == adios2_type_integer8) then
  
            call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)
  
          else
            print*,"===> IO Error: unsuported integer type"
            stop 1
          end if
  
      class default
        print*,"Unknown IO process handler type"
        stop 1

      end select

    end subroutine

    !> @brief Read a 2D 32-bit integer array from file
    !
    !> @todo Check if the "mode" can be read from the configuration file
    !
    !> param[in]    io_proc  : ADIOS2 IO process used for reading
    !> param[in]    var_name : Name of integer array to read
    !> param[in]    start    : What global index to start reading from
    !> param[in]    count    : How many array element to read
    !> param[input] var      : The 2D integer array
    module subroutine read_array_int32_2D(io_proc, var_name, start, count, var)
      class(io_process), intent(in) :: io_proc
      character (len=*), intent(in) :: var_name
      integer(int64), dimension(2), intent(in) :: start
      integer(int64), dimension(2), intent(in) :: count
      integer(int32), dimension(:,:), intent(inout) :: var

      type(adios2_variable):: adios2_var
      integer(int64), dimension(:,:), allocatable :: tmp_var64
      integer(accs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
          call adios2_set_selection(adios2_var, 2, start, count, ierr)

          if (adios2_var%type == adios2_type_integer8) then

            allocate(tmp_var64(size(var, dim=1),size(var,dim=2)))

            call downcast_warning()
            call adios2_get(io_proc%engine, adios2_var, tmp_var64, adios2_mode_sync, ierr)
            var = int(tmp_var64, int32)
  
          else if (adios2_var%type == adios2_type_integer4) then
  
            call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)
  
          else
            print*,"===> IO Error: unsuported integer type"
            stop 1
          end if

      class default
        print*,"Unknown IO process handler type"
        stop 1

      end select

    end subroutine

    !> @brief Read a 2D 64-bit integer array from file
    !
    !> @todo Check if the "mode" can be read from the configuration file
    !
    !> param[in]    io_proc  : ADIOS2 IO process used for reading
    !> param[in]    var_name : Name of integer array to read
    !> param[in]    start    : What global index to start reading from
    !> param[in]    count    : How many array element to read
    !> param[input] var      : The 2D integer array
    module subroutine read_array_int64_2D(io_proc, var_name, start, count, var)
      class(io_process), intent(in) :: io_proc
      character (len=*), intent(in) :: var_name
      integer(int64), dimension(2), intent(in) :: start
      integer(int64), dimension(2), intent(in) :: count
      integer(int64), dimension(:,:), intent(inout) :: var

      type(adios2_variable):: adios2_var
      integer(int32), dimension(:,:), allocatable :: tmp_var32
      integer(accs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
          call adios2_set_selection(adios2_var, 2, start, count, ierr)

          if (adios2_var%type == adios2_type_integer4) then

            allocate(tmp_var32(size(var, dim=1),size(var,dim=2)))

            call adios2_get(io_proc%engine, adios2_var, tmp_var32, adios2_mode_sync, ierr)
            var = int(tmp_var32, int64)
  
          else if (adios2_var%type == adios2_type_integer8) then
  
            call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)
  
          else
            print*,"===> IO Error: unsuported integer type"
            stop 1
          end if

      class default
        print*,"Unknown IO process handler type"
        stop 1

      end select

    end subroutine

    !> @brief Read a 1D 32-bit real array from file
    !
    !> @todo Check if the "mode" can be read from the configuration file
    !
    !> param[in]    io_proc  : ADIOS2 IO process used for reading
    !> param[in]    var_name : Name of real array to read
    !> param[in]    start    : What global index to start reading from
    !> param[in]    count    : How many array element to read
    !> param[input] var      : The 1D real array
    module subroutine read_array_real32_1D(io_proc, var_name, start, count, var)
      class(io_process), intent(in) :: io_proc
      character (len=*), intent(in) :: var_name
      integer(int64), dimension(1), intent(in) :: start
      integer(int64), dimension(1), intent(in) :: count
      real(real32), dimension(:), intent(inout) :: var

      type(adios2_variable):: adios2_var
      real(real64), dimension(:), allocatable :: tmp_var64
      integer(accs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
          call adios2_set_selection(adios2_var, 1, start, count, ierr)

          if (adios2_var%type == adios2_type_dp) then

            allocate(tmp_var64(size(var)))

            call downcast_warning()
            call adios2_get(io_proc%engine, adios2_var, tmp_var64, adios2_mode_sync, ierr)
            var = real(tmp_var64, real32)

          else if (adios2_var%type == adios2_type_real) then

            call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

          else
            print*,"===> IO Error: unsuported real type"
            stop 1
          end if

        class default
          print*,"Unknown IO process handler type"
          stop 1

      end select

    end subroutine

    !> @brief Read a 1D 64-bit real array from file
    !
    !> @todo Check if the "mode" can be read from the configuration file
    !
    !> param[in]    io_proc  : ADIOS2 IO process used for reading
    !> param[in]    var_name : Name of real array to read
    !> param[in]    start    : What global index to start reading from
    !> param[in]    count    : How many array element to read
    !> param[input] var      : The 1D real array
    module subroutine read_array_real64_1D(io_proc, var_name, start, count, var)
      class(io_process), intent(in) :: io_proc
      character (len=*), intent(in) :: var_name
      integer(int64), dimension(1), intent(in) :: start
      integer(int64), dimension(1), intent(in) :: count
      real(real64), dimension(:), intent(inout) :: var

      type(adios2_variable):: adios2_var
      real(real32), dimension(:), allocatable :: tmp_var32
      integer(accs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
          call adios2_set_selection(adios2_var, 1, start, count, ierr)

          if (adios2_var%type == adios2_type_real) then

            allocate(tmp_var32(size(var)))

            call adios2_get(io_proc%engine, adios2_var, tmp_var32, adios2_mode_sync, ierr)
            var = real(tmp_var32, real64)

          else if (adios2_var%type == adios2_type_dp) then

            call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

          else
            print*,"===> IO Error: unsuported real type"
            stop 1
          end if

        class default
          print*,"Unknown IO process handler type"
          stop 1

      end select

    end subroutine

    !> @brief Read a 2D 32-bit real array from file
    !
    !> @todo Check if the "mode" can be read from the configuration file
    !
    !> param[in]    io_proc  : ADIOS2 IO process used for reading
    !> param[in]    var_name : Name of real array to read
    !> param[in]    start    : What global index to start reading from
    !> param[in]    count    : How many array element to read
    !> param[input] var      : The 2D real array
    module subroutine read_array_real32_2D(io_proc, var_name, start, count, var)
      class(io_process), intent(in) :: io_proc
      character (len=*), intent(in) :: var_name
      integer(int64), dimension(2), intent(in) :: start
      integer(int64), dimension(2), intent(in) :: count
      real(real32), dimension(:,:), intent(inout) :: var

      type(adios2_variable):: adios2_var
      real(real64), dimension(:,:), allocatable :: tmp_var64
      integer(accs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)
          call adios2_set_selection(adios2_var, 2, start, count, ierr)

          if (adios2_var%type == adios2_type_dp) then

            allocate(tmp_var64(size(var,dim=1),size(var,dim=2)))

            call downcast_warning()
            call adios2_get(io_proc%engine, adios2_var, tmp_var64, adios2_mode_sync, ierr)
            var = real(tmp_var64, real32)

          else if (adios2_var%type == adios2_type_real) then

            call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

          else
            print*,"===> IO Error: unsuported real type"
            stop 1
          end if

      class default
        print*,"Unknown IO process handler type"
        stop 1

      end select

    end subroutine

    !> @brief Read a 2D 64-bit real array from file
    !
    !> @todo Check if the "mode" can be read from the configuration file
    !
    !> param[in]    io_proc  : ADIOS2 IO process used for reading
    !> param[in]    var_name : Name of real array to read
    !> param[in]    start    : What global index to start reading from
    !> param[in]    count    : How many array element to read
    !> param[input] var      : The 2D real array
        module subroutine read_array_real64_2D(io_proc, var_name, start, count, var)
      class(io_process), intent(in) :: io_proc
      character (len=*), intent(in) :: var_name
      integer(int64), dimension(2), intent(in) :: start
      integer(int64), dimension(2), intent(in) :: count
      real(real64), dimension(:,:), intent(inout) :: var

      type(adios2_variable):: adios2_var
      real(real32), dimension(:,:), allocatable :: tmp_var32
      integer(accs_int) :: ierr

      select type(io_proc)
        type is(adios2_io_process)

          call adios2_inquire_variable(adios2_var, io_proc%io_task, var_name, ierr)


          call adios2_set_selection(adios2_var, 2, start, count, ierr)
          if (adios2_var%type == adios2_type_real) then

            allocate(tmp_var32(size(var,dim=1),size(var,dim=2)))

            call adios2_get(io_proc%engine, adios2_var, tmp_var32, adios2_mode_sync, ierr)
            var = real(tmp_var32, real64)

          else if (adios2_var%type == adios2_type_dp) then

            call adios2_get(io_proc%engine, adios2_var, var, adios2_mode_sync, ierr)

          else
            print*,"===> IO Error: unsuported real type"
            stop 1
          end if

      class default
        print*,"Unknown IO process handler type"
        stop 1
        
      end select

    end subroutine

    !> @brief Print out downcast warning
    subroutine downcast_warning()
      print*,"===> IO Warning:"
      print*,"===> Downcasting from 64-bit to 32-bit, possible loss of precision."
    end subroutine

    end submodule