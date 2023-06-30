!v Module file scalars.mod
!
!  Defines the scalar transport subroutines.

module scalars

  use types, only: field, ccs_mesh, fluid
  use parallel_types, only: parallel_environment

  implicit none

  private

  public :: update_scalars

  interface
    !> Subroutine to perform scalar transport for all scalar fields.
    module subroutine update_scalars(par_env, mesh, flow)
      class(parallel_environment), allocatable, intent(in) :: par_env   !< parallel environment
      type(ccs_mesh), intent(in) :: mesh                                !< the mesh
      type(fluid), intent(inout) :: flow                                !< The structure containting all the fluid fields
    end subroutine update_scalars
  end interface

end module scalars
