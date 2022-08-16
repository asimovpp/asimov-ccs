!v Module file timestepping_mod.f90
!
!  Provides the interface for timestepping methods.
module timestepping

  use kinds, only: ccs_real
  use types, only: field, ccs_mesh, ccs_vector, ccs_matrix

  implicit none

  private

  public :: apply_timestep
  
  !real(ccs_real) :: dt

  interface
    module subroutine apply_timestep(mesh, phi, diag, M, b)
      use mat, only : set_matrix_diagonal
      
      type(ccs_mesh), intent(in) :: mesh
      class(field), intent(in) :: phi
      class(ccs_vector), intent(inout) :: diag
      class(ccs_matrix), intent(inout) :: M
      class(ccs_vector), intent(inout) :: b
    end subroutine
  end interface
end module timestepping
