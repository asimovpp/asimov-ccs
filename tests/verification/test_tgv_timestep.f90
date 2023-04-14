!> @brief Test grid convergence of the 2d tgv test case on cartesian mesh
!
!> @description Generates several square meshes varying dt and get the error to compute their convergence order
program test_tgv_timestep
#include "ccs_macros.inc"

  use testing_lib
  use error_analysis, only: get_order, print_error_summary
  use tgv2d_core, only: run_tgv2d

  implicit none

  type(ccs_mesh), target :: mesh
  integer(ccs_int), parameter :: num_dt = 4
  integer(ccs_int), parameter :: nvar = 3
  real(ccs_real), dimension(nvar, num_dt) :: errors_L2
  real(ccs_real), dimension(nvar, num_dt) :: errors_Linf
  real(ccs_real), dimension(:), allocatable :: orders_L2
  real(ccs_real), dimension(:), allocatable :: orders_Linf
  real(ccs_real), dimension(num_dt) :: refinements
  real(ccs_real), dimension(num_dt) :: dt_list
  real(ccs_real) :: dt

  integer(ccs_int) :: i
  character(len=12), dimension(nvar) :: variable_labels


  call init()

  variable_labels = (/ "U", "V", "P" /)
  dt_list = (/ 0.01, 0.005, 0.0025, 0.00125 /)
  refinements = dt_list / minval(dt_list)

  errors_L2(:, :) = 0.0_ccs_real
  errors_Linf(:, :) = 0.0_ccs_real

  do i = 1, num_dt
    dt = dt_list(i)
    call run_tgv2d(par_env, errors_L2(:, i), errors_Linf(:, i), input_dt=dt)
  end do

  if (par_env%proc_id == par_env%root) then
 
    call print_error_summary(variable_labels, refinements, errors_L2, errors_Linf)

    call get_order(refinements, errors_L2, orders_L2)
    call get_order(refinements, errors_Linf, orders_Linf)

    call assert_gt(orders_L2(1), 0.9_ccs_real, "U not converging in 2nd order ")
    call assert_gt(orders_L2(2), 0.9_ccs_real, "V not converging in 2nd order ")
    !call assert_gt(orders_L2(3), 0.9_ccs_real, "P not converging in 2nd order ")

    call assert_gt(orders_Linf(1), 0.9_ccs_real, "U not converging in 2nd order ")
    call assert_gt(orders_Linf(2), 0.9_ccs_real, "V not converging in 2nd order ")
    !call assert_gt(orders_Linf(3), 0.9_ccs_real, "P not converging in 2nd order ")

  end if

  call fin()

 contains

end program test_tgv_timestep
