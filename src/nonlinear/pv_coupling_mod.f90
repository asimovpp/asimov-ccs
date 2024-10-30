!v Module file pv_coupling.mod
!
!  An interface to pressure-velocity coupling methods (SIMPLE, etc)

module pv_coupling

  use kinds, only: ccs_int, ccs_real
  use types, only: field, ccs_mesh, fluid
  use parallel_types, only: parallel_environment

  implicit none

  private

  public :: solve_nonlinear

  interface

    module subroutine solve_nonlinear(par_env, mesh, eval_sources, it_start, it_end, res_target, &
                                      flow, diverged)
      class(parallel_environment), allocatable, intent(in) :: par_env
      type(ccs_mesh), intent(in) :: mesh
      interface
        !v Subroutine to evaluate source terms, case-specific.
        !
        !  Note this should return the integrated source.
        subroutine eval_sources(flow, phi, R, S)
          use types, only: fluid, field, ccs_vector
          type(fluid), intent(in) :: flow !< Provides access to full flow field
          class(field), intent(in) :: phi !< Field being transported
          class(ccs_vector), intent(inout) :: R !< Work vector (for evaluating linear/implicit sources)
          class(ccs_vector), intent(inout) :: S !< Work vector (for evaluating fixed/explicit sources)
        end subroutine eval_sources
      end interface
      integer(ccs_int), intent(in) :: it_start, it_end
      real(ccs_real), intent(in) :: res_target
      type(fluid), intent(inout) :: flow                              !< Container for flow fields
      logical, optional, intent(out) :: diverged
    end subroutine solve_nonlinear

  end interface

end module pv_coupling
