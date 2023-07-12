!v Submodule file fv_discretisation.smod
!
!  Implementations of the finite volume method using the various discretisation schemes scheme
!
!  @build discretisation

submodule(fv) fv_discretisation

  use vec, only: get_vector_data, restore_vector_data
  use meshing, only: get_face_interpolation, get_local_index
  use types, only: neighbour_locator
  use meshing, only: get_distance, get_centre

  implicit none

contains
  !> Calculates advection coefficient for neighbouring cell using CDS discretisation
  module subroutine calc_advection_coeff_cds(phi, loc_f, mf, bc, coeff)
    type(central_field), intent(in) :: phi !< scalar field
    type(face_locator), intent(in) :: loc_f !< face locator
    real(ccs_real), intent(in) :: mf       !< mass flux at the face
    integer(ccs_int), intent(in) :: bc     !< flag indicating whether cell is on boundary
    real(ccs_real), intent(out) :: coeff   !< advection coefficient to be calculated

    real(ccs_real) :: interpolation_factor

    ! Dummy usage to prevent unused argument.
    associate (scalar => phi)
    end associate
    associate (mflux => mf)
    end associate

    if (bc == 0) then
      call get_face_interpolation(loc_f, interpolation_factor)
      interpolation_factor = 1.0_ccs_real - interpolation_factor
    else
      interpolation_factor = 0.5_ccs_real !1.0_ccs_real
    end if
    coeff = interpolation_factor
  end subroutine calc_advection_coeff_cds

  !> Calculates advection coefficient for neighbouring cell using UDS discretisation
  module subroutine calc_advection_coeff_uds(phi, loc_f, mf, bc, coeff)
    type(upwind_field), intent(in) :: phi !< scalar field
    type(face_locator), intent(in) :: loc_f !< face locator
    real(ccs_real), intent(in) :: mf      !< mass flux at the face
    integer(ccs_int), intent(in) :: bc    !< flag indicating whether cell is on boundary
    real(ccs_real), intent(out) :: coeff  !< advection coefficient to be calculated

    ! Dummy usage to prevent unused argument.
    associate (scalar => phi, foo => bc, bar => loc_f)
    end associate

    if (mf < 0.0) then
      coeff = 1.0_ccs_real
    else
      coeff = 0.0_ccs_real
    end if

  end subroutine calc_advection_coeff_uds

  !> Calculates advection coefficient for neighbouring cell using gamma discretisation
  module subroutine calc_advection_coeff_gamma(phi, loc_f, mf, bc, loc_p, loc_nb, coeff)
    type(gamma_field), intent(inout) :: phi       !< scalar field
    type(face_locator), intent(in) :: loc_f       !< face locator
    real(ccs_real), intent(in) :: mf              !< mass flux at the face
    integer(ccs_int), intent(in) :: bc            !< flag indicating whether cell is on boundary
    type(cell_locator), intent(in) :: loc_p       !< current cell locator
    type(neighbour_locator), intent(in) :: loc_nb !< neighbour cell locator
    real(ccs_real), intent(out) :: coeff          !< advection coefficient to be calculated

    real(ccs_real), dimension(:), pointer :: phi_data
    real(ccs_real), dimension(:), pointer :: dphidx, dphidy, dphidz
    real(ccs_real), dimension(3) :: dphiF, dphiP, d
    real(ccs_real) :: phiF, phiP, dphi, ddphi, phiPt, gamma_m, beta_m

    integer(ccs_int) :: index_p, index_nb, index_bc

    !store values of phi filed in phi_data array
    call get_vector_data(phi%values, phi_data)

    !store x-gradients of phi in dphidx array
    call get_vector_data(phi%x_gradients, dphidx)

    !store y-gradients of phi in dphidx array
    call get_vector_data(phi%y_gradients, dphidy)

    !store z-gradients of phi in dphidx array
    call get_vector_data(phi%z_gradients, dphidz)

    !get the local index of current cell and neighbouring cell
    call get_local_index(loc_p, index_p)
    call get_local_index(loc_nb, index_nb)

    !Dummy usage to prevent unused argument.
    associate (scalar => phi, foo => bc, bar => loc_f)
    end associate

    beta_m = 0.35_ccs_real !value can be varied between 0.1 and 0.5

    if (bc == 0) then
      ! Interior face
      ! -------------
      if (mf < 0.0) then
        phiP = phi_data(index_nb)
        phiF = phi_data(index_p)

        !Gradient of phi at cell center (current cell)
        dphiP(1) = dphidx(index_nb)
        dphiP(2) = dphidy(index_nb)
        dphiP(3) = dphidz(index_nb)

        !Gradient phi at cell face
        dphi = phiF - phiP

        !Get the distance between present and neighbouring cell centers and store it in d
        call get_distance(loc_p, loc_nb, d)
        d = -1.0_ccs_real * d

        !calculate the normalized phi
        ddphi = 2.0_ccs_real * dot_product(dphiP, d)
        phiPt = 1.0_ccs_real - (dphi / ddphi)

        if (phiPt <= 0.0_ccs_real .or. phiPt >= 1.0_ccs_real) then !UD
          coeff = 1.0_ccs_real
        else if (phiPt > beta_m .and. phiPt < 1.0_ccs_real) then !CDS
          coeff = 0.5_ccs_real
        else if (phiPt > 0.0_ccs_real .and. phiPt <= beta_m) then !Gamma
          gamma_m = phiPt / beta_m
          coeff = 1.0_ccs_real - 0.5_ccs_real * gamma_m
        end if
      else
        phiP = phi_data(index_p)
        phiF = phi_data(index_nb)

        !Gradient of phi at cell center (current cell)
        dphiP(1) = dphidx(index_p)
        dphiP(2) = dphidy(index_p)
        dphiP(3) = dphidz(index_p)

        !Gradient phi at cell face
        dphi = phiF - phiP

        !Get the distance between present and neighbouring cell centers and store it in d
        call get_distance(loc_p, loc_nb, d)

        !calculate the normalized phi
        ddphi = 2.0_ccs_real * dot_product(dphiP, d)
        phiPt = 1.0_ccs_real - (dphi / ddphi)

        if (phiPt <= 0.0_ccs_real .or. phiPt >= 1.0_ccs_real) then !UD
          coeff = 0.0_ccs_real
        else if (phiPt > beta_m .and. phiPt < 1.0_ccs_real) then !CDS
          coeff = 0.5_ccs_real
        else if (phiPt > 0.0_ccs_real .and. phiPt <= beta_m) then !Gamma
          gamma_m = phiPt / beta_m
          coeff = 0.5_ccs_real * gamma_m
        end if
      end if

    else
      ! Boundary face
      ! -------------
      ! (Assume CDS - is this correct?)
      coeff = 0.5_ccs_real
    endif

    ! Restore vectors
    call restore_vector_data(phi%values, phi_data)
    call restore_vector_data(phi%x_gradients, dphidx)
    call restore_vector_data(phi%y_gradients, dphidy)
    call restore_vector_data(phi%z_gradients, dphidz)

  end subroutine calc_advection_coeff_gamma
end submodule fv_discretisation
