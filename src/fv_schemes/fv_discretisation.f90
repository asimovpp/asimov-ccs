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
  pure module subroutine calc_advection_coeff_cds(phi, loc_f, mf, bc, coeffaP, coeffaF)
    type(central_field), intent(in) :: phi  !< scalar field
    type(face_locator), intent(in) :: loc_f !< face locator
    real(ccs_real), intent(in) :: mf        !< mass flux at the face
    integer(ccs_int), intent(in) :: bc      !< flag indicating whether cell is on boundary
    real(ccs_real), intent(out) :: coeffaP  !< advection coefficient for current cell
    real(ccs_real), intent(out) :: coeffaF  !< advection coefficient for neighbour cell

    real(ccs_real) :: interpolation_factor

    ! Dummy usage to prevent unused argument.
    associate (scalar => phi)
    end associate
    associate (mflux => mf)
    end associate

    if (bc == 0 .and. (.not. phi%enable_cell_corrections)) then
      call get_face_interpolation(loc_f, interpolation_factor)
      interpolation_factor = 1.0_ccs_real - interpolation_factor
    else
      interpolation_factor = 0.5_ccs_real !1.0_ccs_real
    end if
    coeffaF = interpolation_factor
    coeffaP = 1.0_ccs_real - coeffaF
  end subroutine calc_advection_coeff_cds

  !> Calculates advection coefficient for neighbouring cell using UDS discretisation
  pure module subroutine calc_advection_coeff_uds(phi, loc_f, mf, bc, coeffaP, coeffaF)
    type(upwind_field), intent(in) :: phi   !< scalar field
    type(face_locator), intent(in) :: loc_f !< face locator
    real(ccs_real), intent(in) :: mf        !< mass flux at the face
    integer(ccs_int), intent(in) :: bc      !< flag indicating whether cell is on boundary
    real(ccs_real), intent(out) :: coeffaP  !< advection coefficient for current cell
    real(ccs_real), intent(out) :: coeffaF  !< advection coefficient for neighbour cell

    ! Dummy usage to prevent unused argument.
    associate (scalar => phi, foo => bc, bar => loc_f)
    end associate

    if (mf < 0.0) then
      coeffaF = 1.0_ccs_real
      coeffaP = 1.0_ccs_real - coeffaF
    else
      coeffaF = 0.0_ccs_real
      coeffaP = 1.0_ccs_real - coeffaF
    end if

  end subroutine calc_advection_coeff_uds

  !> Calculates advection coefficient for neighbouring cell using gamma discretisation
  !
  ! The implementation of the Gamma scheme is based on the Dolfyn implementation
  ! https://bazaar.launchpad.net/~hwkrus/dolfyn-cfd/trunk/view/411/src/diffschemes.f90
  ! The original code is distributed under the APACHE-2.0 license reproduced below
  !
  ! Copyright 2003-2014 Henk Krus, Cyclone Fluid Dynamics BV
  ! All Rights Reserved.
  !
  ! Licensed under the Apache License, Version 2.0 (the "License");
  ! you may not use this file except in compliance with the License.
  ! You may obtain a copy of the License at
  !
  ! http://www.dolfyn.net/license.html
  !
  ! Unless required by applicable law or agreed to in writing, software
  ! distributed under the License is distributed on an
  ! "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
  ! either express or implied. See the License for the specific
  ! language governing permissions and limitations under the License.
  !
  module subroutine calc_advection_coeff_gamma(phi, loc_f, mf, bc, loc_p, loc_nb, coeffaP, coeffaF)
    type(gamma_field), intent(inout) :: phi       !< scalar field
    type(face_locator), intent(in) :: loc_f       !< face locator
    real(ccs_real), intent(in) :: mf              !< mass flux at the face
    integer(ccs_int), intent(in) :: bc            !< flag indicating whether cell is on boundary
    type(cell_locator), intent(in) :: loc_p       !< current cell locator
    type(neighbour_locator), intent(in) :: loc_nb !< neighbour cell locator
    real(ccs_real), intent(out) :: coeffaP        !< advection coefficient for current cell
    real(ccs_real), intent(out) :: coeffaF        !< advection coefficient for neighbour cell

    real(ccs_real), dimension(:), pointer :: phi_data
    real(ccs_real), dimension(:), pointer :: dphidx, dphidy, dphidz
    real(ccs_real), dimension(3) :: dphiP, d
    real(ccs_real) :: phiF, phiP, dphi, ddphi, phiPt, gamma_m, beta_m
    real(ccs_real) :: interpolation_factor

    integer(ccs_int) :: index_p, index_nb

    ! store values of phi filed in phi_data array
    call get_vector_data(phi%values, phi_data)

    ! store x-gradients of phi in dphidx array
    call get_vector_data(phi%x_gradients, dphidx)

    ! store y-gradients of phi in dphidx array
    call get_vector_data(phi%y_gradients, dphidy)

    ! store z-gradients of phi in dphidx array
    call get_vector_data(phi%z_gradients, dphidz)

    ! get the local index of current cell and neighbouring cell
    call get_local_index(loc_p, index_p)
    call get_local_index(loc_nb, index_nb)

    ! Dummy usage to prevent unused argument.
    associate (scalar => phi, foo => bc, bar => loc_f)
    end associate

    beta_m = 0.35_ccs_real ! value can be varied between 0.1 and 0.5

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
          coeffaF = 1.0_ccs_real
        else if (phiPt > beta_m .and. phiPt < 1.0_ccs_real) then !CDS
          call get_face_interpolation(loc_f, interpolation_factor)
          interpolation_factor = 1.0_ccs_real - interpolation_factor
          coeffaF = interpolation_factor
        else !Gamma
          gamma_m = phiPt / beta_m
          call get_face_interpolation(loc_f, interpolation_factor)
          coeffaF = 1.0_ccs_real + (gamma_m * (interpolation_factor - 1.0_ccs_real))
        end if
        coeffaP = 1.0_ccs_real - coeffaF
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
          coeffaF = 0.0_ccs_real
        else if (phiPt > beta_m .and. phiPt < 1.0_ccs_real) then !CDS
          call get_face_interpolation(loc_f, interpolation_factor)
          interpolation_factor = 1.0_ccs_real - interpolation_factor
          coeffaF = interpolation_factor
        else !Gamma
          gamma_m = phiPt / beta_m
          call get_face_interpolation(loc_f, interpolation_factor)
          coeffaF = (1.0_ccs_real - interpolation_factor) * gamma_m
        end if
        coeffaP = 1.0_ccs_real - coeffaF
      end if

    else
      ! Boundary face
      ! -------------
      coeffaF = 0.5_ccs_real
      coeffaP = 1.0_ccs_real - coeffaF
    endif

    ! Restore vectors
    call restore_vector_data(phi%values, phi_data)
    call restore_vector_data(phi%x_gradients, dphidx)
    call restore_vector_data(phi%y_gradients, dphidy)
    call restore_vector_data(phi%z_gradients, dphidz)

  end subroutine calc_advection_coeff_gamma

  !> Calculates advection coefficient for neighbouring cell using Linear Upwind discretisation
  !
  ! The implementation of the linear upwind scheme is based on the Dolfyn implementation
  ! https://bazaar.launchpad.net/~hwkrus/dolfyn-cfd/trunk/view/411/src/diffschemes.f90
  ! The original code is distributed under the APACHE-2.0 license reproduced below
  !
  ! Copyright 2003-2014 Henk Krus, Cyclone Fluid Dynamics BV
  ! All Rights Reserved.
  !
  ! Licensed under the Apache License, Version 2.0 (the "License");
  ! you may not use this file except in compliance with the License.
  ! You may obtain a copy of the License at
  !
  ! http://www.dolfyn.net/license.html
  !
  ! Unless required by applicable law or agreed to in writing, software
  ! distributed under the License is distributed on an
  ! "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
  ! either express or implied. See the License for the specific
  ! language governing permissions and limitations under the License.
  !
  module subroutine calc_advection_coeff_luds(phi, loc_f, mf, bc, loc_p, loc_nb, coeffaP, coeffaF)
    type(linear_upwind_field), intent(inout) :: phi     !< scalar field
    type(face_locator), intent(in) :: loc_f       !< face locator
    real(ccs_real), intent(in) :: mf              !< mass flux at the face
    integer(ccs_int), intent(in) :: bc            !< flag indicating whether cell is on boundary
    type(cell_locator), intent(in) :: loc_p       !< current cell locator
    type(neighbour_locator), intent(in) :: loc_nb !< neighbour cell locator
    real(ccs_real), intent(out) :: coeffaP        !< advection coefficient for current cell
    real(ccs_real), intent(out) :: coeffaF        !< advection coefficient for neighbour cell

    real(ccs_real), dimension(:), pointer :: phi_data
    real(ccs_real), dimension(:), pointer :: dphidx, dphidy, dphidz
    real(ccs_real), dimension(3) :: dphiP, d
    real(ccs_real) :: phiF, phiP, dphi, ddphi, phiPt

    integer(ccs_int) :: index_p, index_nb

    ! store values of phi filed in phi_data array
    call get_vector_data(phi%values, phi_data)

    ! store x-gradients of phi in dphidx array
    call get_vector_data(phi%x_gradients, dphidx)

    ! store y-gradients of phi in dphidx array
    call get_vector_data(phi%y_gradients, dphidy)

    ! store z-gradients of phi in dphidx array
    call get_vector_data(phi%z_gradients, dphidz)

    ! get the local index of current cell and neighbouring cell
    call get_local_index(loc_p, index_p)
    call get_local_index(loc_nb, index_nb)

    ! Dummy usage to prevent unused argument.
    associate (scalar => phi, foo => bc, bar => loc_f)
    end associate

    if (bc == 0) then
      ! Internal face
      ! -------------
      if (mf < 0.0) then
        phiP = phi_data(index_nb)
        phiF = phi_data(index_p)

        ! Gradient of phi at cell center (current cell)
        dphiP(1) = dphidx(index_nb)
        dphiP(2) = dphidy(index_nb)
        dphiP(3) = dphidz(index_nb)

        ! Gradient phi at cell face
        dphi = phiF - phiP

        ! Get the distance between present and neighbouring cell centers and store it in d
        call get_distance(loc_p, loc_nb, d)
        d = -1.0_ccs_real * d

        ! calculate the normalized phi
        ddphi = 2.0_ccs_real * dot_product(dphiP, d)
        phiPt = 1.0_ccs_real - (dphi / ddphi)

        if (phiPt <= 0.0_ccs_real .or. phiPt >= 1.0_ccs_real) then ! UD
          coeffaF = 1.0_ccs_real
          coeffaP = 0.0_ccs_real
        else !LUDS
          call get_distance(loc_nb, loc_f, d)
          coeffaF = 1.0_ccs_real
          if (dabs(phiP) > 0.0_ccs_real) then
            coeffaF = coeffaF + (dot_product(dphiP, d) / phiP)
          end if
          coeffaP = 0.0_ccs_real
        end if
      else
        phiP = phi_data(index_p)
        phiF = phi_data(index_nb)

        ! Gradient of phi at cell center (current cell)
        dphiP(1) = dphidx(index_p)
        dphiP(2) = dphidy(index_p)
        dphiP(3) = dphidz(index_p)

        ! Gradient phi at cell face
        dphi = phiF - phiP

        ! Get the distance between present and neighbouring cell centers and store it in d
        call get_distance(loc_p, loc_nb, d)

        ! calculate the normalized phi
        ddphi = 2.0_ccs_real * dot_product(dphiP, d)
        phiPt = 1.0_ccs_real - (dphi / ddphi)

        if (phiPt <= 0.0_ccs_real .or. phiPt >= 1.0_ccs_real) then ! UD
          coeffaF = 0.0_ccs_real
          coeffaP = 1.0_ccs_real
        else !LUDS
          call get_distance(loc_p, loc_f, d)
          coeffaP = 1.0_ccs_real
          if (dabs(phiP) > 0.0_ccs_real) then
            coeffaP = coeffaP + (dot_product(dphiP, d) / phiP)
          end if
          coeffaF = 0.0_ccs_real
        end if
      end if
    else
      ! Boundary face
      ! -------------
      coeffaF = 0.5_ccs_real
      coeffaP = 0.5_ccs_real
    end if

    ! Restore vectors
    call restore_vector_data(phi%values, phi_data)
    call restore_vector_data(phi%x_gradients, dphidx)
    call restore_vector_data(phi%y_gradients, dphidy)
    call restore_vector_data(phi%z_gradients, dphidz)

  end subroutine calc_advection_coeff_luds

end submodule fv_discretisation
