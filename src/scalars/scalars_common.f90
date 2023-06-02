!v Submodule file scalars_common.smod
!
!  Implementation of the scalar transport subroutines

submodule(scalars) scalars_common

contains
  
  !> Subroutine to perform scalar transport for all scalar fields.
  module subroutine calculate_scalars(par_env, mesh, flow)
    class(parallel_environment), allocatable, intent(in) :: par_env   !< parallel environment
    type(ccs_mesh), intent(in) :: mesh                                !< the mesh
    type(fluid), intent(inout) :: flow                                !< The structure containting all the fluid fields

    ! XXX: Temporarily silence unused variables
    associate(foo => par_env, bar => mesh, baz => flow)
    end associate
    
  end subroutine calculate_scalars
  
end submodule scalars_common
