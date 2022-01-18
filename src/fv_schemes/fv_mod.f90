!> @brief Module file fv.mod
!
!> @details An interface to finite volume implementations (CDS, UDS, etc.)

module fv

  use kinds, only : accs_real, accs_int
  use types, only : matrix, vector, mesh

  implicit none

  private

  public :: compute_fluxes
  public :: calc_advection_coeff
  public :: calc_mass_flux
  public :: calc_cell_coords

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

	! Calculates advection coefficient for neighbouring cell 
    module subroutine calc_advection_coeff(ngb_idx, self_idx, face_area, coeff, cps, u, v, BC)
      integer(accs_int), intent(in) :: ngb_idx, self_idx
      real(accs_real), intent(in) :: face_area
      real(accs_real), intent(inout) :: coeff
      integer(accs_int), intent(in) :: cps
      real(accs_real), dimension(:,:) :: u, v
      integer(accs_int), intent(in) :: BC
    end subroutine calc_advection_coeff

	! Calculates mass flux across given edge. Note: assuming rho = 1 and uniform grid
    module function calc_mass_flux(edge_len, u, v, ngb_row, ngb_col, self_row, self_col, BC_flag) result(flux)
      real(accs_real), intent(in) :: edge_len
      real(accs_real), dimension(:,:), intent(in) :: u, v
      integer(accs_int), intent(in) :: ngb_row, ngb_col
      integer(accs_int), intent(in) :: self_row, self_col
      integer(accs_int), intent(in) :: BC_flag
      real(accs_real) :: flux
    end function calc_mass_flux

    ! Assigns source vector
    ! Calculates the row and column indices from flattened vector index
    ! Note: assumes square mesh
    module subroutine calc_cell_coords(idx, cps, row, col)
      integer(accs_int), intent(in) :: idx, cps
      integer(accs_int), intent(out) :: row, col
    end subroutine calc_cell_coords

  end interface

end module fv