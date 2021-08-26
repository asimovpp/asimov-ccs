module compute

  use iso_fortran_env
  use parallel_types, only: parallel_environment

  implicit none

  private

  interface

  module subroutine compute_pi(num_steps, par_env, mypi)
    integer(kind=int64) , intent(in) :: num_steps
    class(parallel_environment), intent(in) :: par_env
    double precision, intent(out) :: mypi
  end subroutine

  end interface 

  public :: compute_pi

end module compute