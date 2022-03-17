!> @brief A simple test to see if the basic functionality of init and fin is in place.
program init_fin
  use testing_lib
  implicit none

  call init()
  call fin()
end program init_fin
