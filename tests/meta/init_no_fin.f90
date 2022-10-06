!v Init should fail if fin has not been called.
program init_no_fin
  use testing_lib
  implicit none

  call init()
end program init_no_fin
