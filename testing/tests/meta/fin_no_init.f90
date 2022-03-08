!> @brief Fin should fail if init has not been called.
program fin_no_init
  use testing_lib
  implicit none

  call fin()
end program fin_no_init
