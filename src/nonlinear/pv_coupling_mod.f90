!v Module file pv_coupling.mod
!
!  An interface to pressure-velocity coupling methods (SIMPLE, etc)

module pv_coupling

  use core, only: ccs_options
  use types, only: field, ccs_mesh, fluid
  use parallel_types, only: parallel_environment

  implicit none

  private

  public :: solve_nonlinear

  interface

    module subroutine solve_nonlinear(par_env, run_options, mesh, flow, diverged)
      class(parallel_environment), allocatable, intent(in) :: par_env
      type(ccs_options), intent(in) :: run_options
      type(ccs_mesh), intent(in) :: mesh
      type(fluid), intent(inout) :: flow                              !< Container for flow fields
      logical, optional, intent(out) :: diverged
    end subroutine solve_nonlinear

  end interface

end module pv_coupling
