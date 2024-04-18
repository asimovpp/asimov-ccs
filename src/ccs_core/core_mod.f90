!v Module file core.mod
!
! The core module implements the core functionality to define and run a CCS case.

module core

  use types, only: ccs_options
  use parallel_types, only: parallel_environment

  implicit none

  private
  public :: initialise_mesh

  interface

    module subroutine initialise_mesh(par_env, shared_env, run_options)
      class(parallel_environment), intent(in) :: par_env
      class(parallel_environment), intent(in) :: shared_env
      type(ccs_options), intent(in) :: run_options
    end subroutine initialise_mesh

  end interface
end module core
