!v Program file for 2D TaylorGreenVortex case
!
!  @build mpi+petsc

program tgv2d
#include "ccs_macros.inc"

  use tgv2d_core
  use constants, only: ccs_split_type_low_high
  use parallel, only: initialise_parallel_environment, create_new_par_env 

  implicit none

  class(parallel_environment), allocatable :: par_env
  class(parallel_environment), allocatable :: shared_env
  real(ccs_real), dimension(4) :: error_L2
  real(ccs_real), dimension(4) :: error_Linf
  logical :: use_mpi_splitting

  ! Launch MPI
  call initialise_parallel_environment(par_env)
  use_mpi_splitting = .false.
  call create_new_par_env(par_env, ccs_split_type_low_high, use_mpi_splitting, shared_env)

  call run_tgv2d(par_env, shared_env, error_L2, error_Linf)

  ! Finalise MPI
  call cleanup_parallel_environment(par_env)

end program tgv2d
