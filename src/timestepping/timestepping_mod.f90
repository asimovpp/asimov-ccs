!v Module file timestepping_mod.f90
!
!  Provides the interface for timestepping methods.
module timestepping

  use kinds, only: ccs_real
  use types, only: field, ccs_mesh, ccs_vector, ccs_matrix

  implicit none

  private

  public :: apply_timestep
  public :: set_timestep
  public :: get_timestep
  public :: update_old_values
  
  interface
    module subroutine apply_timestep(mesh, phi, diag, M, b)
      type(ccs_mesh), intent(in) :: mesh
      class(field), intent(in) :: phi
      class(ccs_vector), intent(inout) :: diag
      class(ccs_matrix), intent(inout) :: M
      class(ccs_vector), intent(inout) :: b
    end subroutine

    module subroutine set_timestep(timestep)
      real(ccs_real), intent(in) :: timestep
    end subroutine
    
    module function get_timestep() result(timestep)
      real(ccs_real) :: timestep
    end function 

    module subroutine update_old_values(x)
      class(field), intent(inout) :: x 
    end subroutine

  end interface

end module timestepping
