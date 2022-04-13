!> @brief Module file fv.mod
!
!> @details An interface to finite volume implementations (CDS, UDS, etc.)

module fv

  use kinds, only : ccs_real, ccs_int
  use types, only : ccs_matrix, ccs_vector, ccs_mesh, field, upwind_field, central_field, bc_config, face_locator

  implicit none

  private

  public :: compute_fluxes
  public :: calc_advection_coeff
  public :: calc_diffusion_coeff
  public :: calc_mass_flux
  public :: calc_cell_coords
  public :: update_gradient
  
  interface calc_advection_coeff
    module procedure calc_advection_coeff_cds
    module procedure calc_advection_coeff_uds
  end interface calc_advection_coeff

  interface

  !> @brief Calculates advection coefficient for neighbouring cell using CDS discretisation
  !
  !> @param[in] phi         - scalar (central) field
  !> @param[in] mf          - mass flux at the face
  !> @param[in] bc          - flag indicating whether cell is on boundary
  !> @param[out] coeff      - advection coefficient to be calculated
  module subroutine calc_advection_coeff_cds(phi, mf, bc, coeff)
    type(central_field), intent(in) :: phi
    real(ccs_real), intent(in) :: mf
    integer(ccs_int), intent(in) :: bc
    real(ccs_real), intent(out) :: coeff
  end subroutine calc_advection_coeff_cds
  
  !> @brief Calculates advection coefficient for neighbouring cell using UDS discretisation
  !
  !> @param[in] phi         - scalar (upwind) field
  !> @param[in] mf          - mass flux at the face
  !> @param[in] bc          - flag indicating whether cell is on boundary
  !> @param[out] coeff      - advection coefficient to be calculated
  module subroutine calc_advection_coeff_uds(phi, mf, bc, coeff)
    type(upwind_field), intent(in) :: phi
    real(ccs_real), intent(in) :: mf
    integer(ccs_int), intent(in) :: bc
    real(ccs_real), intent(out) :: coeff
  end subroutine calc_advection_coeff_uds

  !> @brief Sets the diffusion coefficient
  !
  !> @param[in] index_p - the local cell index
  !> @param[in] index_nb  - the local neigbouring cell index
  !> @param[in] mesh      - the mesh structure
  !> @param[out] coeff         - the diffusion coefficient
  !
  ! XXX: why is this a function when the equivalent advection ones are subroutines?
  module function calc_diffusion_coeff(index_p, index_nb, mesh) result(coeff)
    integer(ccs_int), intent(in) :: index_p
    integer(ccs_int), intent(in) :: index_nb
    type(ccs_mesh), intent(in) :: mesh
    real(ccs_real) :: coeff
  end function calc_diffusion_coeff

  !> @brief Computes fluxes and assign to matrix and RHS
  !
  !> @param[in] phi       - scalar field structure
  !> @param[in] mf        - mass flux field structure (defined at faces)
  !> @param[in] mesh - the mesh being used
  !> @param[in] bcs       - the boundary conditions structure being used
  !> @param[in] cps       - the number of cells per side in the (square) mesh
  !> @param[in,out] M     - Data structure containing matrix to be filled
  !> @param[in,out] vec   - Data structure containing RHS vector to be filled
  module subroutine compute_fluxes(phi, mf, mesh, bcs, cps, M, vec)
    class(field), intent(in) :: phi
    class(field), intent(in) :: mf
    type(ccs_mesh), intent(in) :: mesh
    type(bc_config), intent(in) :: bcs
    integer(ccs_int), intent(in) :: cps
    class(ccs_matrix), intent(inout) :: M
    class(ccs_vector), intent(inout) :: vec   
  end subroutine


  !> @brief Calculates mass flux across given face. Note: assumes rho = 1 and uniform grid
  !
  !> @param[in] u, v     - arrays containing x, y velocities
  !> @param[in] p        - array containing pressure
  !> @param[in] dpdx   - array containing pressure gradient in x
  !> @param[in] dpdy   - array containing pressure gradient in y
  !> @param[in] invAx    - array containing inverse momentum diagonal in x
  !> @param[in] invAy    - array containing inverse momentum diagonal in y
  !> @param[in] loc_f    - face locator
  !> @param[out] flux    - The flux across the boundary
  module function calc_mass_flux(u, v, p, dpdx, dpdy, invAu, invAv, loc_f) result(flux)
    real(ccs_real), dimension(:), intent(in) :: u, v
    real(ccs_real), dimension(:), intent(in) :: p
    real(ccs_real), dimension(:), intent(in) :: dpdx, dpdy
    real(ccs_real), dimension(:), intent(in) :: invAu, invAv
    type(face_locator), intent(in) :: loc_f
    real(ccs_real) :: flux
  end function calc_mass_flux

  !> @brief Calculates the row and column indices from flattened vector index. Assumes square mesh
  !
  !> @param[in] index  - cell index
  !> @param[in] cps    - number of cells per side
  !> @param[out] row   - cell row within mesh
  !> @param[out] col   - cell column within mesh
  module subroutine calc_cell_coords(index, cps, row, col)
    integer(ccs_int), intent(in) :: index, cps
    integer(ccs_int), intent(out) :: row, col
  end subroutine calc_cell_coords

  !> @brief Performs an update of the gradients of a field.
  !
  !> @param[in]    mesh - the mesh
  !> @param[inout] phi       - the field whose gradients we want to update
  !
  !> @note This will perform a parallel update of the gradient fields to ensure halo cells are
  !!       correctly updated on other PEs.
  module subroutine update_gradient(mesh, phi)
    type(ccs_mesh), intent(in) :: mesh
    class(field), intent(inout) :: phi
  end subroutine update_gradient
  
  end interface

end module fv
