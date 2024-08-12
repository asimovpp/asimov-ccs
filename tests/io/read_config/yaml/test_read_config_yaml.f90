!v Test program for the yaml parser.

program test_read_config_yaml

  use testing_lib

  use fortran_yaml_c_interface, only: parse
  use read_config, only: get_value
  
  implicit none

  character(len=*), parameter :: config_yaml = "config.yml"
  class(*), pointer :: config_file !< Pointer to config file
  
  logical :: expected, required
  
  call init()
  call setup()

  expected = .true. ! Expect to find values
  required = .true. ! Require to find values
  call test_read_values(config_file, expected, required)
  call test_read_dict(config_file)
  
  call fin()

contains

  ! Perform test-case setup
  subroutine setup()

    character(len=:), allocatable :: error

    ! Open config file
    config_file => parse(config_yaml, error)
    if (allocated(error)) then
      call stop_test(trim(error))
    end if

  end subroutine setup

  subroutine test_read_values(conf_file, expected, required)

    class(*), pointer, intent(in) :: conf_file
    logical, intent(in) :: expected !< Is the value expected to be found?
    logical, intent(in) :: required !< Is the value required to be found?

    integer :: ival
    real(ccs_real) :: rval
    character(len=:), allocatable :: sval
    logical :: lval

    logical :: val_present

    print *, "- read integer"
    call get_value(conf_file, "integer", ival, val_present, required=required)
    if (val_present .neqv. expected) then
      call stop_test("Integer value status does not meet expectation")
    end if
    if (val_present) then
      if (ival /= 42) then
        call stop_test("Integer value does not match expectation")
      end if
    end if

    print *, "- read real"
    call get_value(conf_file, "real", rval, val_present, required=required)
    if (val_present .neqv. expected) then
      call stop_test("Real value status does not meet expectation")
    end if
    if (val_present) then
      if (rval /= 3.14_ccs_real) then
        call stop_test("Real value does not match expectation")
      end if
    end if

    print *, "- read string"
    call get_value(conf_file, "string", sval, val_present, required=required)
    if (val_present .neqv. expected) then
      call stop_test("String value status does not meet expectation")
    end if
    if (val_present) then
      if (sval /= "Hello world!") then
        call stop_test("String value does not match expectation")
      end if
    end if

    print *, "- read logical true"
    call get_value(conf_file, "yes", lval, val_present, required=required)
    if (val_present .neqv. expected) then
      call stop_test("Logical (true) value status does not meet expectation")
    end if
    if (val_present) then
      if (.not. lval) then
        call stop_test("Logial (true) value does not match expectation")
      end if
    end if
    print *, "- read logical false"
    call get_value(conf_file, "no", lval, val_present, required=required)
    if (val_present .neqv. expected) then
      call stop_test("Logical (false) value status does not meet expectation")
    end if
    if (val_present) then
      if (lval) then
        call stop_test("Logial (false) value does not match expectation")
      end if
    end if
    
  end subroutine test_read_values

  ! Test reading values from a dictionary
  subroutine test_read_dict(conf_file)

    use fortran_yaml_c, only: type_dictionary, type_error

    class(*), pointer, intent(in) :: conf_file

    class(*), pointer :: dict
    type(type_error), allocatable :: io_err

    logical :: expected, required

    character(len=:), allocatable :: sval
    logical :: val_present

    select type(conf_file)
    type is(type_dictionary)

      print *, "Reading dictionary values"

      ! Read from dictionary with values
      print *, "+ Dictionary: present"
      dict => conf_file%get_dictionary("present", required=.true., error=io_err)
      expected = .true. ! Expect to find values
      required = .true. ! Require to find values
      call test_read_values(dict, expected, required)

      !- Read an optional string
      print *, "- read optional string"
      call get_value(dict, "opt", sval, val_present, required=.false.)
      if (.not. val_present) then
        call stop_test("Optional string value status does not meet expectation")
      end if
      if (val_present) then
        if (sval /= "HI!") then
          print *, sval
          call stop_test("Optional string value does not match expectation")
        end if
      end if

      ! Read from dictionary without values
      print *, "+ Dictionary: absent"
      dict => conf_file%get_dictionary("absent", required=.true., error=io_err)
      expected = .false. ! Expect to find values
      required = .false. ! Require to find values
      call test_read_values(dict, expected, required)

      ! Try to read a value as a dictionary
      dict => conf_file%get_dictionary("not_a_dict", required=.true., error=io_err)
      if (.not. allocated(io_err)) then
        call stop_test("Non-dictionary should have raised an error")
      end if

      ! Read from optional dictionary
      print *, "+ Optional dictionary"
      dict => conf_file%get_dictionary("absent", required=.false., error=io_err)
      call get_value(dict, "foo", sval, val_present, required=.false.)
      if (.not. val_present) then
        call stop_test("Optional dictionary string value status does not meet expectation")
      end if
      if (val_present) then
        if (sval /= "bar") then
          call stop_test("Optional dictionary string value does not match expectation")
        end if
      end if

    class default
      call stop_test("Invalid config file dictionary")
    end select

  end subroutine test_read_dict
  
end program test_read_config_yaml
