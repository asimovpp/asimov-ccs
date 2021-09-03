module compute

  use iso_fortran_env
  use parallel_types, only: parallel_environment

  implicit none

  private

  public :: compute_pi

  interface

  !> @brief Compute Pi
  !>
  !> @param[in] integer(kind=int64) num_steps - number of steps for the algorithm
  !> @param[in] parallel_environment par_env - the parallel environment
  !> @param[out] double precision mypi - computed value for Pi
  module subroutine compute_pi(num_steps, par_env, mypi)
    integer(kind=int64) , intent(in) :: num_steps
    class(parallel_environment), intent(in) :: par_env
    double precision, intent(out) :: mypi
  end subroutine

  end interface 

end module compute
