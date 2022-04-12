!> @brief A simple test to see if the basic functionality of stop_test is in place.
program call_stop_test
  use testing_lib
  implicit none

  call init()

  call stop_test("Stopping test.")

  call fin()
end program call_stop_test
