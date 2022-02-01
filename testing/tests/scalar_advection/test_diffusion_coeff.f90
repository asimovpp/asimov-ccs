!> @brief Test that cells have correct numbers of neighbours
!
!> @description for any mesh with >1 cell, every cell must have at least 1 neighbour.
program test_mesh_neighbours

  use testing_lib
  use mesh_utils, only : build_square_mesh
  use fv, only: calc_diffusion_coeff

  type(mesh) :: square_mesh
  integer(accs_int), parameter :: cps = 50
  real(accs_real) :: coeff
  real(accs_real), parameter :: expected_coeff = 2.e-2

  
  call init()

  square_mesh = build_square_mesh(cps, 1.0_accs_real, par_env)

  coeff = calc_diffusion_coeff(1,1,square_mesh)

  if (abs(coeff - expected_coeff) .le. tiny(coeff)) then
    write(message, *) "FAIL: incorrect diffusion coefficient computed. Expected ", expected_coeff, " computed ", coeff
    call stop_test(message)
  end if
  
  call fin()

end program test_mesh_neighbours
