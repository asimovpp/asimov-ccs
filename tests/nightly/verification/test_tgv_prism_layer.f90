!> @brief Test grid convergence of the 2d tgv test case that includes a prism layer
!
!> @description Generates several square meshes varying cps and get the error to compute their convergence order
program test_tgv_prism_layer
#include "ccs_macros.inc"

  use testing_lib
  use error_analysis, only: get_order, print_error_summary
  use mesh_utils, only: build_square_mesh
  use tgv2d_core, only: run_tgv2d, domain_size
  use mesh_utils, only: compute_face_interpolation
  use meshing, only: get_local_num_cells

  implicit none

  type(ccs_mesh), target :: mesh
  integer(ccs_int), parameter :: num_cps = 5
  integer(ccs_int), parameter :: nvar = 3
  real(ccs_real), dimension(nvar, num_cps) :: error_L2
  real(ccs_real), dimension(nvar, num_cps) :: error_Linf
  real(ccs_real), dimension(:), allocatable :: orders_L2
  real(ccs_real), dimension(:), allocatable :: orders_Linf
  real(ccs_real), dimension(:), allocatable :: min_error_L2
  real(ccs_real), dimension(:), allocatable :: min_error_Linf
  real(ccs_real), dimension(num_cps) :: refinements
  real(ccs_real) :: growth_rate
  integer(ccs_int), dimension(num_cps) :: cps_list
  integer(ccs_int) :: cps

  character(len=12), dimension(nvar) :: variable_labels

  integer(ccs_int) :: i, j

  call init()

  variable_labels = (/"U", "V", "P"/)
  domain_size = 3.14159265358979323
  growth_rate = 1.2

  cps_list = (/16, 32, 64, 128, 256/)
  refinements = real(maxval(cps_list(:))) / real(cps_list(:))

  error_L2(:, :) = 0.0_ccs_real
  error_Linf(:, :) = 0.0_ccs_real

  do i = 1, num_cps
    cps = cps_list(i)
    mesh = build_square_mesh(par_env, cps, domain_size)

    call generate_prism_layer(growth_rate, cps, mesh)

    call run_tgv2d(par_env, error_L2(:, i), error_Linf(:, i), mesh)
  end do

  if (par_env%proc_id == par_env%root) then

    call print_error_summary(variable_labels, refinements, error_L2, error_Linf)

    call get_order(refinements, error_L2, orders_L2)
    call get_order(refinements, error_Linf, orders_Linf)

    call assert_gt(orders_L2(1), 1.9_ccs_real, "U not converging in 2nd order ")
    call assert_gt(orders_L2(2), 1.9_ccs_real, "V not converging in 2nd order ")
    !call assert_gt(orders_L2(3), 1.9_ccs_real, "P not converging in 2nd order ")

    call assert_gt(orders_Linf(1), 1.4_ccs_real, "U not converging in 2nd order ")
    call assert_gt(orders_Linf(2), 1.4_ccs_real, "V not converging in 2nd order ")
    !call assert_gt(orders_Linf(3), 1.4_ccs_real, "P not converging in 2nd order ")

  end if

  call fin()

contains

  !v Modifies a square mesh by applying a growth rate along the x axis hence generating a prism layer
  subroutine generate_prism_layer(growth_rate, cps, mesh)
    real(ccs_real), intent(in) :: growth_rate
    integer(ccs_int), intent(in) :: cps
    type(ccs_mesh), intent(inout) :: mesh

    real(ccs_real) :: power, dx
    integer(ccs_int) :: iface, icell, local_num_cells

    dx = domain_size / real(cps)

    ! Update face centers
    mesh%geo%x_f(1, :, :) = apply_gr(mesh%geo%x_f(1, :, :), growth_rate)

    ! Update vertices
    mesh%geo%vert_coords(1, :, :) = apply_gr(mesh%geo%vert_coords(1, :, :), growth_rate)

    ! Update face areas
    call get_local_num_cells(mesh, local_num_cells)
    do icell = 1, local_num_cells
      do iface = 1, 4
        if (abs(mesh%geo%face_normals(2, iface, icell)) .gt. 0.01) then
          mesh%geo%face_areas(iface, icell) = apply_gr(mesh%geo%x_p(1, icell) + dx/2, growth_rate) - apply_gr(mesh%geo%x_p(1, icell) - dx/2, growth_rate)
        end if
      end do
    end do

    ! Update cell volumes
    mesh%geo%volumes(:) = dx * (apply_gr(mesh%geo%x_p(1, :) + dx/2, growth_rate) - apply_gr(mesh%geo%x_p(1, :) - dx/2, growth_rate))

    ! Update cell centers (average of the new faces positions)
    mesh%geo%x_p(1, :) = 0.5_ccs_real * (apply_gr(mesh%geo%x_p(1, :) + dx/2, growth_rate) + apply_gr(mesh%geo%x_p(1, :) - dx/2, growth_rate))

    ! Update face interpolation
    call compute_face_interpolation(mesh)

  end subroutine

  !v applies a growth rate to the x coordinate
  elemental function apply_gr(x, growth_rate) result(new_x)
    real(ccs_real), intent(in) :: x
    real(ccs_real), intent(in) :: growth_rate
    real(ccs_real) :: new_x
    real(ccs_real) :: power

    power = log(growth_rate + 1.0_ccs_real) / log(2.0_ccs_real)
    new_x = domain_size * (x / domain_size)**power

  end function

end program test_tgv_prism_layer
