submodule(fv_kernels:advection_common) linear_upwind_advection
  use types
  use vec, only: get_vector_data, restore_vector_data
  use meshing, only: get_face_interpolation, get_local_index, get_distance, get_centre

  implicit none

contains

  module procedure advect_luw_coeffs
  real(ccs_real) :: coeffaP        !< advection coefficient for current cell
  real(ccs_real) :: coeffaF        !< advection coefficient for neighbour cell

  real(ccs_real), dimension(:), pointer :: phi_data
  real(ccs_real), dimension(:), pointer :: dphidx, dphidy, dphidz
  real(ccs_real), dimension(3) :: dphiP, d
  real(ccs_real) :: phiF, phiP, dphi, ddphi, phiPt

  integer(ccs_int) :: index_p, index_nb

  ! store values of phi filed in phi_data array
  call get_vector_data(phi % values, phi_data)

  ! store x-gradients of phi in dphidx array
  call get_vector_data(phi % x_gradients, dphidx)

  ! store y-gradients of phi in dphidx array
  call get_vector_data(phi % y_gradients, dphidy)

  ! store z-gradients of phi in dphidx array
  call get_vector_data(phi % z_gradients, dphidz)

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

  coeffs = [coeffaP, coeffaF]

  ! Restore vectors
  call restore_vector_data(phi % values, phi_data)
  call restore_vector_data(phi % x_gradients, dphidx)
  call restore_vector_data(phi % y_gradients, dphidy)
  call restore_vector_data(phi % z_gradients, dphidz)
  end procedure advect_luw_coeffs

end submodule advect_luw_impl
