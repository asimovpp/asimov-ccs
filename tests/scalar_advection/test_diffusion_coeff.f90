!> @brief Test that the diffusion coefficient is being calculated correctly
!
!> @description Currently hard-coded result, waiting for better treatment of diffusion
program test_diffusion_coeff

  use testing_lib
  use mesh_utils, only : build_square_mesh
  use fv, only: calc_diffusion_coeff

  type(mesh) :: square_mesh
  integer(ccs_int), parameter :: cps = 50
  real(ccs_real) :: coeff
  real(ccs_real), parameter :: expected_coeff = -1.e-2_ccs_real

  
  call init()

  square_mesh = build_square_mesh(par_env, cps, 1.0_ccs_real)

  coeff = calc_diffusion_coeff(1,1,square_mesh)

  if (abs(coeff - expected_coeff) .ge. tiny(coeff)) then
    write(message, *) "FAIL: incorrect diffusion coefficient computed. Expected ", expected_coeff, " computed ", coeff
    call stop_test(message)
  end if
  
  call fin()

end program test_diffusion_coeff
