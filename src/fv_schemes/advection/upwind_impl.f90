submodule(fv_kernels:advection_mod) upwind_advection_impl
implicit none

contains
!> Calculates advection coefficient for neighbouring cell using UDS discretisation
module procedure advect_upwind_coeffs
real(ccs_real) :: coeffaP  !< advection coefficient for current cell
real(ccs_real) :: coeffaF  !< advection coefficient for neighbour cell

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
