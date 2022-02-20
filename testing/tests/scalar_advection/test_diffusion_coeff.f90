!> @brief Test that the diffusion coefficient is being calculated correctly
!
!> @description Currently hard-coded result, waiting for better treatment of diffusion
program test_diffusion_coeff

  use testing_lib
  use mesh_utils, only : build_square_mesh
  use fv, only: calc_diffusion_coeff

  type(mesh) :: square_mesh
  integer(accs_int), parameter :: cps = 50
  real(accs_real) :: coeff
  real(accs_real), parameter :: expected_coeff = -2.e-2_accs_real

  
  call init()

  square_mesh = build_square_mesh(cps, 1.0_accs_real, par_env)

  coeff = calc_diffusion_coeff(1,1,square_mesh)

  if (abs(coeff - expected_coeff) .ge. tiny(coeff)) then
    write(message, *) "FAIL: incorrect diffusion coefficient computed. Expected ", expected_coeff, " computed ", coeff
    call stop_test(message)
  end if
  
  call fin()

end program test_diffusion_coeff
