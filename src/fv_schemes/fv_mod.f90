!> @brief Module file fv.mod
!
!> @details An interface to finite volume implementations (CDS, UDS, etc.)

module fv

  use kinds, only : accs_real, accs_int
  use types, only : matrix, vector, mesh, field, upwind_field, central_field

  implicit none

  private

  public :: compute_fluxes
  public :: calc_advection_coeff
  public :: calc_mass_flux
  public :: calc_cell_coords

  !interface calc_advection_coeff
  !  module procedure calc_advection_coeff_cds
  !  module procedure calc_advection_coeff_uds
  !end interface calc_advection_coeff

  interface

  module subroutine calc_advection_coeff(ngb_idx, self_idx, face_area, coeff, cps, u, v, BC)
    integer(accs_int), intent(in) :: ngb_idx, self_idx
    real(accs_real), intent(in) :: face_area
    real(accs_real), intent(inout) :: coeff
    integer(accs_int), intent(in) :: cps
    class(field), intent(in) :: u, v
    integer(accs_int), intent(in) :: BC
  end subroutine calc_advection_coeff
    
  !> @brief Calculates advection coefficient for neighbouring cell using the CDS scheme
  !
  !> @param[in] ngb_idx - neighbour index
  !> @param[in] self_idx - cell index
  !> @param[in] face_area - area of face between cell and neighbour
  !> @param[in,out] coeff - advection coefficient that is computed
  !> @param[in] cps - number of cells per side
  !> @param[in] u, v - velocity fields in x, y directions
  !> @param[in] BC - flag to indicate boundary
  module subroutine calc_advection_coeff_cds(ngb_idx, self_idx, face_area, coeff, cps, u, v, BC)
    integer(accs_int), intent(in) :: ngb_idx, self_idx
    real(accs_real), intent(in) :: face_area
    real(accs_real), intent(inout) :: coeff
    integer(accs_int), intent(in) :: cps
    type(central_field), intent(in) :: u, v
    integer(accs_int), intent(in) :: BC
  end subroutine calc_advection_coeff_cds
  
  !> @brief Calculates advection coefficient for neighbouring cell using the UDS scheme
  !
  !> @param[in] ngb_idx - neighbour index
  !> @param[in] self_idx - cell index
  !> @param[in] face_area - area of face between cell and neighbour
  !> @param[in,out] coeff - advection coefficient that is computed
  !> @param[in] cps - number of cells per side
  !> @param[in] u, v - velocity fields in x, y directions
  !> @param[in] BC - flag to indicate boundary
  module subroutine calc_advection_coeff_uds(ngb_idx, self_idx, face_area, coeff, cps, u, v, BC)
    integer(accs_int), intent(in) :: ngb_idx, self_idx
    real(accs_real), intent(in) :: face_area
    real(accs_real), intent(inout) :: coeff
    integer(accs_int), intent(in) :: cps
    type(upwind_field), intent(in) :: u, v
    integer(accs_int), intent(in) :: BC
  end subroutine calc_advection_coeff_uds

    !> @brief Computes fluxes and assign to matrix and RHS
    !
    !> @param[in,out] mat - Data structure containing matrix to be filled
    !> @param[in,out] vec - Data structure containing RHS vector to be filled
    !> @param[in] u, v - arrays containing velocity fields in x, y directions
    !> @param[in] cell_mesh - the mesh being used
    module subroutine compute_fluxes(mat, vec, u, v, cell_mesh)
      class(matrix), intent(inout) :: mat
      class(vector), intent(inout) :: vec   
      class(field), intent(in) :: u, v
      type(mesh), intent(in) :: cell_mesh
    end subroutine


    !> @brief Calculates mass flux across given edge. Note: assumes rho = 1 and uniform grid
    !
    !> @param[in] edge_len - length of edge
    !> @param[in] u, v - velocity field in x, y directions
    !> @param[in] ngb_row, ngb_col - row and column index of neighbouring cell
    !> @param[in] self_row, self_col - row and column index of cell
    !> @param[in] BC_flag - indicates whether a cell is on a boundary and which boundary it is.
    module function calc_mass_flux(edge_len, u, v, ngb_row, ngb_col, self_row, self_col, BC_flag) result(flux)
      real(accs_real), intent(in) :: edge_len
      class(field), intent(in) :: u, v
      integer(accs_int), intent(in) :: ngb_row, ngb_col
      integer(accs_int), intent(in) :: self_row, self_col
      integer(accs_int), intent(in) :: BC_flag
      real(accs_real) :: flux
    end function calc_mass_flux

    !> @brief Calculates the row and column indices from flattened vector index
    !
    !> @param[in] idx - cell index
    !> @param[in] cps - number of cells per side in mesh
    !> @param[out] row, col - the row and column of the cell in the mesh
    module subroutine calc_cell_coords(idx, cps, row, col)
      integer(accs_int), intent(in) :: idx, cps
      integer(accs_int), intent(out) :: row, col
    end subroutine calc_cell_coords

  end interface

end module fv
