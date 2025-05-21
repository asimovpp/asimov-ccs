submodule(fv_kernels:advection_common) central_difference_advection
  use types
  use vec, only: get_vector_data, restore_vector_data
  use meshing, only: get_face_interpolation, get_local_index, get_distance, get_centre

  implicit none

contains
  !> Calculates advection coefficient for neighbouring cell using CDS discretisation
  module procedure advect_cds_coeffs
  real(ccs_real) :: coeffaF
  real(ccs_real) :: coeffaP
  real(ccs_real) :: interpolation_factor

  ! Dummy usage to prevent unused argument.
  associate (scalar => phi)
  end associate
  associate (mflux => mf)
  end associate

  if (bc == 0 .and. (.not. phi % enable_cell_corrections)) then
    call get_face_interpolation(loc_f, interpolation_factor)
    interpolation_factor = 1.0_ccs_real - interpolation_factor
  else
    interpolation_factor = 0.5_ccs_real !1.0_ccs_real
  end if
  coeffaF = interpolation_factor
  coeffaP = 1.0_ccs_real - coeffaF

  coeffs = [coeffaP, coeffaF]
  end procedure advect_cds_coeffs

end submodule advection_cds_impl
