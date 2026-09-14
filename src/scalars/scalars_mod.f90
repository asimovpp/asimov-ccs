!v Module file scalars.mod
!
!  Defines the scalar transport subroutines.

module scalars

  use types, only: field, ccs_mesh, fluid, ccs_residuals
  use parallel_types, only: parallel_environment
  use kinds, only: ccs_real

  implicit none

  private

  public :: update_scalars

  interface
    !> Subroutine to perform scalar transport for all scalar fields.
    module subroutine update_scalars(par_env, mesh, eval_sources, flow, residuals)
      class(parallel_environment), allocatable, intent(in) :: par_env   !< parallel environment
      type(ccs_mesh), intent(in) :: mesh                                !< the mesh
      interface
        subroutine eval_sources(flow, phi, R, S)
          use types, only: fluid, field, ccs_vector
          type(fluid), intent(in) :: flow !< Provides access to full flow field
          class(field), intent(in) :: phi !< Field being transported
          class(ccs_vector), intent(inout) :: R !< Work vector (for evaluating linear/implicit sources)
          class(ccs_vector), intent(inout) :: S !< Work vector (for evaluating fixed/explicit sources)
        end subroutine eval_sources
      end interface
      type(fluid), intent(inout) :: flow                                !< The structure containting all the fluid fields
      type(ccs_residuals), optional, intent(inout) :: residuals
    end subroutine update_scalars
  end interface

end module scalars
