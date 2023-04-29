!v Program file for 2D TaylorGreenVortex case
!
!  @build mpi+petsc

program tgv2d
#include "ccs_macros.inc"

  use tgv2d_core

  implicit none

  class(parallel_environment), allocatable :: par_env
  real(ccs_real), dimension(4) :: error_L2
  real(ccs_real), dimension(4) :: error_Linf

  ! Launch MPI
  call initialise_parallel_environment(par_env)

  call run_tgv2d(par_env, error_L2, error_Linf)

  ! Finalise MPI
  call cleanup_parallel_environment(par_env)

end program tgv2d
