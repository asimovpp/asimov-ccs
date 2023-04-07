!> @brief Test the face interpolation with a 2 cells mesh
!
!> @description Generates a basic mesh geoemtry and topology for a 2 cells mesh and makes sure the face interpolation is 
!!              properly computed and accessed (with get_face_interpolation)
program test_tgv_cartesian
#include "ccs_macros.inc"

  use testing_lib
  use mesh_utils, only: build_square_mesh
  use tgv2d_core, only: run_tgv2d, domain_size

  implicit none

  type(ccs_mesh), target :: mesh
  integer(ccs_int), parameter :: num_cps = 4
  real(ccs_real), dimension(4, num_cps) :: error_L2
  real(ccs_real), dimension(4, num_cps) :: error_Linf
  integer(ccs_int), dimension(num_cps) :: cps_list
  integer(ccs_int) :: cps

  integer(ccs_int) :: i


  call init()

  domain_size = 3.14159265358979323

  cps_list = (/ 4, 8, 16, 32 /)
  error_L2(:, :) = 0.0_ccs_real
  error_Linf(:, :) = 0.0_ccs_real

  do i = 1, num_cps
    cps = cps_list(i)
    mesh = build_square_mesh(par_env, cps, domain_size)

    call run_tgv2d(par_env, error_L2(:, i), error_Linf(:, i), mesh)

  end do

  print *, error_L2
  print *, error_Linf
  !call analyse_error(cps_list, error_L2, slope_L2, min_error_L2)
  !call analyse_error(cps_list, error_Linf, slope_Linf, min_error_Linf)

  !call assert_eq(1.0_ccs_real, 2.0_ccs_real, "test")


  call fin()

! contains

  !subroutine analyse_error(cps_list, errors, slope, min)


  !end subroutine

end program test_tgv_cartesian
