!> @brief Test grid convergence of the 2d tgv test case on cartesian mesh
!
!> @description Generates several square meshes varying cps and get the error to compute their convergence order
program test_tgv_loop
#include "ccs_macros.inc"

  use testing_lib
  use error_analysis, only: get_order, print_error_summary
  use mesh_utils, only: build_square_mesh
  use tgv2d_core, only: run_tgv2d, domain_size

  implicit none

  integer(ccs_int), parameter :: num_cps = 3
  integer(ccs_int), parameter :: nvar = 3
  real(ccs_real), dimension(nvar, num_cps) :: error_L2
  real(ccs_real), dimension(nvar, num_cps) :: error_Linf
  real(ccs_real), dimension(num_cps) :: refinements
  integer(ccs_int), dimension(num_cps) :: cps_list
  integer(ccs_int) :: cps

  character(len=12), dimension(nvar) :: variable_labels

  integer(ccs_int) :: i

  call init()

  variable_labels = (/"U", "V", "P"/)
  domain_size = 3.14159265358979323

  cps_list = (/4, 6, 8/)
  refinements = real(maxval(cps_list(:))) / real(cps_list(:))

  error_L2(:, :) = 0.0_ccs_real
  error_Linf(:, :) = 0.0_ccs_real

  do i = 1, num_cps
    cps = cps_list(i)
    mesh = build_square_mesh(par_env, shared_env, cps, domain_size)

    call run_tgv2d(par_env, shared_env, error_L2(:, i), error_Linf(:, i), input_mesh=mesh)
  end do

  if (par_env%proc_id == par_env%root) then

    call print_error_summary(variable_labels, refinements, error_L2, error_Linf)

  end if

  call fin()

end program test_tgv_loop
