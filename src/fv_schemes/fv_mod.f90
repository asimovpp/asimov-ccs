!> @brief Module file fv.mod
!
!> @details An interface to finite volume implementations (CDS, UDS, etc.)

module fv

  use constants, only : ndim
  use kinds, only : accs_real, accs_int
  use types, only : matrix, vector, mesh, field, upwind_field, central_field, bc_config
  use bc_constants

  implicit none

  private

  public :: compute_fluxes
  public :: calc_advection_coeff
  public :: calc_diffusion_coeff
  public :: calc_mass_flux
  public :: calc_cell_coords

  interface calc_advection_coeff
    module procedure calc_advection_coeff_cds
    module procedure calc_advection_coeff_uds
  end interface calc_advection_coeff

  interface

  !> @brief Calculates advection coefficient for neighbouring cell using CDS discretisation
  !
  !> @param[in] phi         - scalar field
  !> @param[in] mf          - mass flux at the face
  !> @param[in] bc          - flag indicating whether cell is on boundary
  !> @param[out] coeff      - advection coefficient to be calculated
  module subroutine calc_advection_coeff_cds(phi, mf, bc, coeff)
    type(central_field), intent(in) :: phi
    real(accs_real), intent(in) :: mf
    integer(accs_int), intent(in) :: bc
    real(accs_real), intent(out) :: coeff
  end subroutine calc_advection_coeff_cds
  
  !> @brief Calculates advection coefficient for neighbouring cell using UDS discretisation
  !
  !> @param[in] phi         - scalar field
  !> @param[in] mf          - mass flux at the face
  !> @param[in] bc          - flag indicating whether cell is on boundary
  !> @param[out] coeff      - advection coefficient to be calculated
  module subroutine calc_advection_coeff_uds(phi, mf, bc, coeff)
    type(upwind_field), intent(in) :: phi
    real(accs_real), intent(in) :: mf
    integer(accs_int), intent(in) :: bc
    real(accs_real), intent(out) :: coeff
  end subroutine calc_advection_coeff_uds

  !> @brief Sets the diffusion coefficient
  !
  !> @param[in] local_self_idx - the local cell index
  !> @param[in] local_ngb_idx  - the local neigbouring cell index
  !> @param[in] cell_mesh      - the mesh structure
  !> @param[out] coeff         - the diffusion coefficient
  !
  ! XXX: why is this a function when the equivalent advection ones are subroutines?
  module function calc_diffusion_coeff(local_self_idx, local_ngb_idx, cell_mesh) result(coeff)
    integer(accs_int), intent(in) :: local_self_idx
    integer(accs_int), intent(in) :: local_ngb_idx
    type(mesh), intent(in) :: cell_mesh
    real(accs_real) :: coeff
  end function calc_diffusion_coeff

  !> @brief Computes fluxes and assign to matrix and RHS
  !
  !> @param[in] phi       - scalar field structure
  !> @param[in] u, v      - field structures in x, y directions
  !> @param[in] cell_mesh - the mesh being used
  !> @param[in] bcs       - the boundary conditions structure being used
  !> @param[in] cps       - the number of cells per side in the (square) mesh
  !> @param[in,out] M     - Data structure containing matrix to be filled
  !> @param[in,out] vec   - Data structure containing RHS vector to be filled
  module subroutine compute_fluxes(phi, u, v, cell_mesh, bcs, cps, M, vec)
    class(field), intent(in) :: phi
    class(field), intent(in) :: u, v
    type(mesh), intent(in) :: cell_mesh
    type(bc_config), intent(in) :: bcs
    integer(accs_int), intent(in) :: cps
    class(matrix), allocatable, intent(inout) :: M
    class(vector), intent(inout) :: vec   
  end subroutine


  !> @brief Calculates mass flux across given face. Note: assumes rho = 1 and uniform grid
  !
  !> @param[in] u, v     - arrays containing x, y velocities
  !> @param[in] ngb_idx  - Row and column of given neighbour in mesh
  !> @param[in] self_idx - Row and column of given cell in mesh
  !> @param[in] bc_flag  - Flag to indicate if neighbour is a boundary cell
  !> @param[out] flux    - The flux across the boundary
  module function calc_mass_flux(u, v, ngb_idx, self_idx, face_normal, bc_flag) result(flux)
    real(accs_real), dimension(:), intent(in) :: u, v
    integer(accs_int), intent(in) :: ngb_idx
    integer(accs_int), intent(in) :: self_idx
    real(accs_real), dimension(ndim), intent(in) :: face_normal
    integer(accs_int), intent(in) :: bc_flag
    real(accs_real) :: flux
  end function calc_mass_flux

  !> @brief Calculates the row and column indices from flattened vector index. Assumes square mesh
  !
  !> @param[in] idx  - cell index
  !> @param[in] cps  - number of cells per side
  !> @param[out] row - cell row within mesh
  !> @param[out] col - cell column within mesh
  module subroutine calc_cell_coords(idx, cps, row, col)
    integer(accs_int), intent(in) :: idx, cps
    integer(accs_int), intent(out) :: row, col
  end subroutine calc_cell_coords

  end interface

end module fv
