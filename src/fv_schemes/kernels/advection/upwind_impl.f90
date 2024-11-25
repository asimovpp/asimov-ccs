submodule(advection_upwind_mod) advection_upwind_impl
  use advection_upwind_mod
  use vec, only: get_vector_data, restore_vector_data
  use meshing, only: get_face_interpolation, get_local_index, get_distance, get_centre
  use types, only: neighbour_locator

  implicit none

contains
  !> Calculates advection coefficient for neighbouring cell using UDS discretisation
  module procedure advect_upwind_coeffs
  real(ccs_real) :: coeffaP  !< advection coefficient for current cell
  real(ccs_real) :: coeffaF  !< advection coefficient for neighbour cell

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

  coeffs = [coeffaF, coeffaP]

  end procedure advect_upwind_coeffs
end submodule advection_upwind_impl
