!> @brief Test grid convergence of the 2d tgv test case on cartesian mesh
!
!> @description Generates several square meshes varying cps and get the error to compute their convergence order
program test_tgv_cartesian
#include "ccs_macros.inc"

  use testing_lib
  use mesh_utils, only: build_square_mesh
  use tgv2d_core, only: run_tgv2d, domain_size

  implicit none

  type(ccs_mesh), target :: mesh
  integer(ccs_int), parameter :: num_cps = 4
  integer(ccs_int), parameter :: nvar = 3
  real(ccs_real), dimension(nvar, num_cps) :: error_L2
  real(ccs_real), dimension(nvar, num_cps) :: error_Linf
  real(ccs_real), dimension(:), allocatable :: orders_L2
  real(ccs_real), dimension(:), allocatable :: orders_Linf
  real(ccs_real), dimension(:), allocatable :: min_error_L2
  real(ccs_real), dimension(:), allocatable :: min_error_Linf
  integer(ccs_int), dimension(num_cps) :: cps_list
  integer(ccs_int) :: cps

  integer(ccs_int) :: i, j


  call init()

  domain_size = 3.14159265358979323

  cps_list = (/ 8, 16, 32, 64 /)
  error_L2(:, :) = 0.0_ccs_real
  error_Linf(:, :) = 0.0_ccs_real

  do i = 1, num_cps
    cps = cps_list(i)
    mesh = build_square_mesh(par_env, cps, domain_size)

    call run_tgv2d(par_env, error_L2(:, i), error_Linf(:, i), mesh)
  end do

  print *, "-------------------------------------"
  print *, "Summary of errors"

  do j=1, nvar
    if (j == 1) print *, "---- U errors"
    if (j == 2) print *, "---- V errors"
    if (j == 3) print *, "---- P errors"
    do i=1, num_cps
      print *, cps_list(i), error_L2(j, i), error_Linf(j, i)
    end do
  end do

  call analyse_error(cps_list, error_L2, orders_L2, min_error_L2)
  call analyse_error(cps_list, error_Linf, orders_Linf, min_error_Linf)

  print *, "-------------------------------------"
  print *, "Results"

  print *, "U order: ", orders_L2(1), orders_Linf(1)
  print *, "V order: ", orders_L2(2), orders_Linf(2)
  print *, "P order: ", orders_L2(3), orders_Linf(3)

  call assert_gt(orders_L2(1), 1.9_ccs_real, "U not converging in 2nd order ")
  call assert_gt(orders_L2(2), 1.9_ccs_real, "V not converging in 2nd order ")
  call assert_gt(orders_L2(3), 1.9_ccs_real, "P not converging in 2nd order ")

  call assert_gt(orders_Linf(1), 1.9_ccs_real, "U not converging in 2nd order ")
  call assert_gt(orders_Linf(2), 1.9_ccs_real, "V not converging in 2nd order ")
  call assert_gt(orders_Linf(3), 1.9_ccs_real, "P not converging in 2nd order ")

  call fin()

 contains

  subroutine analyse_error(cps_list, errors, orders, min_error)

    integer(ccs_int), dimension(:), intent(in) :: cps_list
    real(ccs_real), dimension(:, :), intent(in) :: errors
    real(ccs_real), dimension(:), allocatable, intent(out) :: orders
    real(ccs_real), dimension(:), allocatable, intent(out) :: min_error
    real(ccs_real), dimension(num_cps) :: x, y
    real(ccs_real) :: x_bar, y_bar, Sxx, Sxy, alpha, beta
    integer(ccs_int) :: i, j, nvar

    nvar = size(errors, dim=1)
    allocate(orders(nvar))
    allocate(min_error(nvar))
    orders(:) = 0.0_ccs_real
    min_error(:) = 0.0_ccs_real

    do i=1, nvar

      x(:) = log(real(maxval(cps_list(:)))/real(cps_list(:)))
      
      y(:) = log(errors(i, :))

      x_bar = sum(x(:))/num_cps
      y_bar = sum(y(:))/num_cps

      Sxy = 0.0_ccs_real
      Sxx = 0.0_ccs_real
      do j=1, num_cps
        Sxy = Sxy + (x(j) - x_bar) * (y(j) - y_bar)
        Sxx = Sxx + (x(j) - x_bar)**2
      end do

      beta = Sxy/Sxx
      alpha = y_bar - beta*x_bar

      orders(i) = beta
      min_error(i) = minval(errors(i, :)) !alpha

    end do

  end subroutine

end program test_tgv_cartesian
