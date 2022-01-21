!> @brief Submodule file fv_CDS.smod
!> @build CDS
!
!> @details An implementation of the finite volume method using the CDS scheme

submodule (fv) fv_CDS

  implicit none

contains
  
  ! Calculates advection coefficient for neighbouring cell 
  module subroutine calc_advection_coeff_cds(ngb_idx, self_idx, face_area, coeff, cps, u, v, BC)
    integer(accs_int), intent(in) :: ngb_idx, self_idx
    real(accs_real), intent(in) :: face_area
    real(accs_real), intent(inout) :: coeff
    integer(accs_int), intent(in) :: cps
    class(central_field), intent(in) :: u, v
    integer(accs_int), intent(in) :: BC

    integer(accs_int) :: ngb_row, ngb_col       ! neighbour coordinates within grid
    integer(accs_int) :: self_row, self_col     ! cell coordinates within grid

    ! Find where we are in the grid first
    call calc_cell_coords(ngb_idx, cps, ngb_row, ngb_col)
    call calc_cell_coords(self_idx, cps, self_row, self_col)

    coeff = calc_mass_flux(face_area, u, v, ngb_row, ngb_col, self_row, self_col, BC)
  end subroutine calc_advection_coeff_cds
  
  ! Calculates advection coefficient for neighbouring cell 
  module subroutine calc_advection_coeff_uds(ngb_idx, self_idx, face_area, coeff, cps, u, v, BC)
    integer(accs_int), intent(in) :: ngb_idx, self_idx
    real(accs_real), intent(in) :: face_area
    real(accs_real), intent(inout) :: coeff
    integer(accs_int), intent(in) :: cps
    class(upwind_field), intent(in) :: u, v
    integer(accs_int), intent(in) :: BC

    integer(accs_int) :: ngb_row, ngb_col       ! neighbour coordinates within grid
    integer(accs_int) :: self_row, self_col     ! cell coordinates within grid

    ! Find where we are in the grid first
    call calc_cell_coords(ngb_idx, cps, ngb_row, ngb_col)
    call calc_cell_coords(self_idx, cps, self_row, self_col)

    coeff = calc_mass_flux(face_area, u, v, ngb_row, ngb_col, self_row, self_col, BC)
    coeff = min(coeff, 0.0_accs_real)
  end subroutine calc_advection_coeff_uds

end submodule fv_CDS
