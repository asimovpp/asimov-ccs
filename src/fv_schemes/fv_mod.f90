!> @brief Module file fv.mod
!
!> @details An interface to finite volume implementations (CDS, UDS, etc.)

module fv

  use kinds, only : accs_real, accs_int
  use types, only : matrix, vector, mesh

  implicit none

  private

  public :: compute_fluxes

  interface
    
	!> @brief Interface to compute fluxes and assign to matrix and RHS
	!
	!> @param[in,out] mat - Data structure containing matrix to be filled
	!> @param[in,out] vec - Data structure containing RHS vector to be filled
	!> @param[in] u, v - arrays containing velocity fields in x, y directions
	!> @param[in] cell_mesh - the mesh being used
	module subroutine compute_fluxes(mat, vec, u, v, cell_mesh)
	  class(matrix), intent(inout) :: mat
      class(vector), intent(inout) :: vec   
      real(accs_real), dimension(:,:), intent(in) :: u, v
	  type(mesh), intent(in) :: cell_mesh
	end subroutine

  end interface

end module fv