!> @brief Test grid convergence of the 2d tgv test case that includes a prism layer
!
!> @description Generates several square meshes varying cps and get the error to compute their convergence order
program test_tgv_prism_layer
#include "ccs_macros.inc"

  use testing_lib
  use ccs_base, only: bnd_names_default
  use error_analysis, only: get_order, print_error_summary
  use types, only: cell_locator, face_locator, vert_locator
  use mesh_utils, only: build_square_mesh
  use tgv2d_core, only: run_tgv2d, domain_size
  use mesh_utils, only: compute_face_interpolation
  use meshing, only: get_local_num_cells, get_total_num_cells, set_mesh_object, nullify_mesh_object, &
                     create_cell_locator, create_face_locator, create_vert_locator, &
                     get_centre, set_centre, set_area, set_volume, get_face_normal, get_volume, get_face_area
  use utils, only: debug_print, str

  implicit none

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
    mesh = build_square_mesh(par_env, shared_env, cps, domain_size, &
         bnd_names_default(1:4))

    call generate_prism_layer(shared_env, growth_rate, cps, mesh)

    call run_tgv2d(par_env, shared_env, error_L2(:, i), error_Linf(:, i), input_mesh=mesh)
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
  subroutine generate_prism_layer(shared_env, growth_rate, cps, mesh)
    class(parallel_environment), intent(in) :: shared_env
    real(ccs_real), intent(in) :: growth_rate
    integer(ccs_int), intent(in) :: cps
    type(ccs_mesh), intent(inout) :: mesh

    real(ccs_real) :: power, dx
    integer(ccs_int) :: iface, icell, ivert, local_num_cells, total_num_cells
    real(ccs_real), dimension(3) :: x_p, x_f, x_v, normal
    real(ccs_real) :: area, volume
    type(cell_locator) :: loc_p
    type(face_locator) :: loc_f
    type(vert_locator) :: loc_v

    call set_mesh_object(mesh)
    
    dx = domain_size / real(cps)
    
    call get_local_num_cells(local_num_cells)
    call get_total_num_cells(total_num_cells)

    ! Update face centers
    do icell = 1, local_num_cells
      call create_cell_locator(icell, loc_p)
      do iface = 1, 4
        call create_face_locator(icell, iface, loc_f)
        call get_centre(loc_f, x_f)
        x_f(1) = apply_gr(x_f(1), growth_rate)
        call set_centre(loc_f, x_f)
      end do

      ! Update vertices
      do ivert = 1, 4
        call create_vert_locator(icell, ivert, loc_v)
        call get_centre(loc_v, x_v)
        x_v(1) = apply_gr(x_v(1), growth_rate)
        call set_centre(loc_v, x_v)
      end do

      ! Update face areas
      do iface = 1, 4
        call create_face_locator(icell, iface, loc_f)
        call get_face_normal(loc_f, normal)
        if (abs(normal(2)) .gt. 0.01) then
          call get_centre(loc_p, x_p)
          area = apply_gr(x_p(1) + dx/2, growth_rate) - apply_gr(x_p(1) - dx/2, growth_rate)
          call set_area(area, loc_f)
        end if
      end do
    end do

    ! Update cell volumes
    do icell = 1, total_num_cells
      call create_cell_locator(icell, loc_p)
      call get_centre(loc_p, x_p)
      volume = dx * (apply_gr(x_p(1) + dx/2, growth_rate) - apply_gr(x_p(1) - dx/2, growth_rate))
      call set_volume(volume, loc_p)

      ! Update cell centers (average of the new faces positions)
      x_p(1) = 0.5_ccs_real * (apply_gr(x_p(1) + dx/2, growth_rate) + apply_gr(x_p(1) - dx/2, growth_rate))
      call set_centre(loc_p, x_p)
    end do

    call sync(shared_env)

    ! Update face interpolation
    call compute_face_interpolation(mesh)

    call nullify_mesh_object()
    
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
