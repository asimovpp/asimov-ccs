!> @brief Submodule file fv_discretisation.smod
!> @build discretisation
!
!> @details Implementations of the finite volume method using the various discretisation schemes scheme

submodule (fv) fv_discretisation

  implicit none

contains
  
  !> @brief Calculates advection coefficient for neighbouring cell using CDS discretisation
  !
  !> @param[in] ngb_idx   - global cell neighbour index
  !> @param[in] self_idx  - global cell index
  !> @param[in] face_area - area of face between self and neighbour
  !> @param[in,out] coeff - advection coefficient to be calculated
  !> @param[in] cps       - number of cells per side in (square) mesh
  !> @param[in] u, v      - velocity fields in x, y directions
  !> @param[in] bc        - flag indicating whether cell is on boundary
  module subroutine calc_advection_coeff_cds(ngb_idx, self_idx, face_area, cps, u, v, bc, coeff)
    integer(accs_int), intent(in) :: ngb_idx, self_idx
    real(accs_real), intent(in) :: face_area
    integer(accs_int), intent(in) :: cps
    type(central_field), intent(in) :: u, v
    integer(accs_int), intent(in) :: bc
    real(accs_real), intent(out) :: coeff

    integer(accs_int) :: ngb_row, ngb_col       ! neighbour coordinates within grid
    integer(accs_int) :: self_row, self_col     ! cell coordinates within grid
    real(accs_real) :: interpolation_factor

    ! Find where we are in the grid first
    call calc_cell_coords(ngb_idx, cps, ngb_row, ngb_col)
    call calc_cell_coords(self_idx, cps, self_row, self_col)

    if (bc == 0) then
      interpolation_factor = 0.5_accs_real
    else
      interpolation_factor = 1.0_accs_real
    end if
    coeff = calc_mass_flux(u, v, ngb_row, ngb_col, self_row, self_col, face_area, bc) * interpolation_factor
  end subroutine calc_advection_coeff_cds
  
  !> @brief Calculates advection coefficient for neighbouring cell using UDS discretisation
  !
  !> @param[in] ngb_idx   - global cell neighbour index
  !> @param[in] self_idx  - global cell index
  !> @param[in] face_area - area of face between self and neighbour
  !> @param[in,out] coeff - advection coefficient to be calculated
  !> @param[in] cps       - number of cells per side in (square) mesh
  !> @param[in] u, v      - velocity fields in x, y directions
  !> @param[in] bc        - flag indicating whether cell is on boundary
  module subroutine calc_advection_coeff_uds(ngb_idx, self_idx, face_area, cps, u, v, bc, coeff)
    integer(accs_int), intent(in) :: ngb_idx, self_idx
    real(accs_real), intent(in) :: face_area
    integer(accs_int), intent(in) :: cps
    type(upwind_field), intent(in) :: u, v
    integer(accs_int), intent(in) :: bc
    real(accs_real), intent(out) :: coeff

    integer(accs_int) :: ngb_row, ngb_col       ! neighbour coordinates within grid
    integer(accs_int) :: self_row, self_col     ! cell coordinates within grid

    ! Find where we are in the grid first
    call calc_cell_coords(ngb_idx, cps, ngb_row, ngb_col)
    call calc_cell_coords(self_idx, cps, self_row, self_col)

    coeff = calc_mass_flux(u, v, ngb_row, ngb_col, self_row, self_col, face_area, bc)
    coeff = min(coeff, 0.0_accs_real)
  end subroutine calc_advection_coeff_uds

end submodule fv_discretisation
