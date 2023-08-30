!v Program file for Poiseuille case
!
!  @build mpi+petsc

program poiseuille
#include "ccs_macros.inc"

  use Poiseuille_core

  implicit none

  class(parallel_environment), allocatable :: par_env
  real(ccs_real), dimension(3) :: error_L2
  real(ccs_real), dimension(3) :: error_Linf

  ! Launch MPI
  call initialise_parallel_environment(par_env)

  call run_poiseuille(par_env, error_L2, error_Linf)

  ! Finalise MPI
  call cleanup_parallel_environment(par_env)

end program poiseuille

