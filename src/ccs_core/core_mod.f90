!v Module file core.mod
!
! The core module implements the core functionality to define and run a CCS case.

module core

  implicit none

  private
  public: initialise_mesh

  interface

    module subroutine initialise_mesh(par_env, shared_env, run_options)
      type(parallel_environment), intent(in) :: par_env
      type(parallel_environment), intent(in) :: shared_env
      type(ccs_options), intent(in) :: run_options
    end subroutine initialise_mesh

  end interface
end module core
