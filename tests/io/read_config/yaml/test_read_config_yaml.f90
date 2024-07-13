!v Test program for the yaml parser.

program test_read_config_yaml

  use testing_lib

  use fortran_yaml_c_interface, only: parse
  use read_config, only: get_value
  
  implicit none

  character(len=*), parameter :: config_yaml = "config.yml"
  class(*), pointer :: config_file !< Pointer to config file
  
  call init()
  call setup()

  call test_read_values(config_file)
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

  subroutine test_read_values(conf_file)

    class(*), pointer, intent(in) :: conf_file

    integer :: ival
    real(ccs_real) :: rval
    character(len=:), allocatable :: sval
    logical :: lval

    logical :: val_present

    call get_value(conf_file, "integer", ival, val_present, required=.true.)
    if (.not. val_present) then
      call stop_test("Integer value was not found")
    end if
    if (ival /= 42) then
      call stop_test("Integer value does not match expectation")
    end if

    call get_value(conf_file, "real", rval, val_present, required=.true.)
    if (.not. val_present) then
      call stop_test("Real value was not found")
    end if
    if (rval /= 3.14_ccs_real) then
      call stop_test("Real value does not match expectation")
    end if

    call get_value(conf_file, "string", sval, val_present, required=.true.)
    if (.not. val_present) then
      call stop_test("String value was not found")
    end if
    if (sval /= "Hello world!") then
      call stop_test("String value does not match expectation")
    end if
    
    call get_value(conf_file, "yes", lval, val_present, required=.true.)
    if (.not. val_present) then
      call stop_test("Logical (true) value was not found")
    end if
    if (.not. lval) then
      call stop_test("Logial (true) value does not match expectation")
    end if
    call get_value(conf_file, "no", lval, val_present, required=.true.)
    if (.not. val_present) then
      call stop_test("Logical (false) value was not found")
    end if
    if (lval) then
      call stop_test("Logial (false) value does not match expectation")
    end if
    
  end subroutine test_read_values

  ! Test reading values from a dictionary
  subroutine test_read_dict(conf_file)

    use fortran_yaml_c, only: type_dictionary, type_error

    class(*), pointer, intent(in) :: conf_file

    class(*), pointer :: dict
    type(type_error), allocatable :: io_err

    select type(conf_file)
    type is(type_dictionary)

      print *, "Reading dictionary values"
      dict => conf_file%get_dictionary("present", required=.true., error=io_err)

      call test_read_values(dict)

    class default
      call stop_test("Invalid config file dictionary")
    end select

  end subroutine test_read_dict
  
end program test_read_config_yaml
