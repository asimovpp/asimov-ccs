!v Set of tools used to analyse discretisation errors
module error_analysis
#include "ccs_macros.inc"

  use kinds
  use types

contains

  !v Prints errors provided in argument in a human readable way, plus displays associated orders
  subroutine print_error_summary(variable_labels, refinements, errors, errors_secondary)
    character(len=12), dimension(:), intent(in) :: variable_labels !< label associated to each error
    real(ccs_real), dimension(:), intent(in) :: refinements !< refinement (likely in time or space) against which the orders are computed
    real(ccs_real), dimension(:, :), intent(in) :: errors !< error values, error(variable, refinement)
    real(ccs_real), dimension(:, :), intent(in), optional :: errors_secondary !< An optional second set of errors to display (could be a different norm for example)

    real(ccs_real), dimension(:), allocatable :: orders
    real(ccs_real), dimension(:), allocatable :: orders_secondary
    character(len=30) :: fmt
    integer(ccs_int) :: nref, nvar, i, j

    nvar = size(errors, dim=1)
    nref = size(refinements)

    call get_order(refinements, errors, orders)
    if (present(errors_secondary)) then
      call get_order(refinements, errors_secondary, orders_secondary)
    end if

    print *, "----------------------------------------------------"
    print *, "Summary of errors"

    do j = 1, nvar
      print *, "------ " // trim(variable_labels(j)) // " errors"
      do i = 1, nref
        if (present(errors_secondary)) then
          fmt = '(f12.4,e12.4,e12.4)'
          write (*, fmt) refinements(i), errors(j, i), errors_secondary(j, i)
        else
          fmt = '(f12.4,e12.4)'
          write (*, fmt) refinements(i), errors(j, i)
        end if
      end do
    end do

    print *, "----------------------------------------------------"
    print *, "Convergence orders"
    do j = 1, nvar
      if (present(errors_secondary)) then
        fmt = '(a12,f12.4,f12.4)'
        write (*, fmt) trim(variable_labels(j)) // " order: ", orders(j), orders_secondary(j)
      else
        fmt = '(a12,f12.4)'
        write (*, fmt) trim(variable_labels(j)) // " order: ", orders(j)
      end if
    end do
    print *, "----------------------------------------------------"
    print *, ""

  end subroutine

  !v Computes convergence orders from a refinement list and the associated errors
  ! The orders are computed as the slope of the linear regression of the log of error against the log of refinements
  pure subroutine get_order(refinements, errors, orders)

    real(ccs_real), dimension(:), intent(in) :: refinements !< refinement (likely in time or space) against which the orders are computed
    real(ccs_real), dimension(:, :), intent(in) :: errors !< error values, error(variable, refinement)
    real(ccs_real), dimension(:), allocatable, intent(out) :: orders !< the computed orders for each variable
    real(ccs_real), dimension(:), allocatable :: x, y
    real(ccs_real) :: x_bar, y_bar, Sxx, Sxy, alpha, beta
    integer(ccs_int) :: nref, nvar, i, j

    nvar = size(errors, dim=1)
    nref = size(refinements)

    allocate (orders(nvar))
    allocate (x(nref))
    allocate (y(nref))

    orders(:) = 0.0_ccs_real

    do i = 1, nvar

      x(:) = log(refinements(:))
      y(:) = log(errors(i, :))

      x_bar = sum(x(:)) / nref
      y_bar = sum(y(:)) / nref

      Sxy = 0.0_ccs_real
      Sxx = 0.0_ccs_real
      do j = 1, nref
        Sxy = Sxy + (x(j) - x_bar) * (y(j) - y_bar)
        Sxx = Sxx + (x(j) - x_bar)**2
      end do

      beta = Sxy / Sxx
      alpha = y_bar - beta * x_bar

      orders(i) = beta

    end do

  end subroutine

  !v Modifies a square mesh by slightly moving cell centres
  subroutine disturb_cartesian(cps, domain_size, mesh)
    use mesh_utils, only: compute_face_interpolation
    use meshing, only: get_total_num_cells

    integer(ccs_int), intent(in) :: cps
    real(ccs_real), intent(in) :: domain_size
    type(ccs_mesh), intent(inout) :: mesh

    real(ccs_real) :: dx
    integer(ccs_int) :: icell, total_num_cells, idim, icell_global
    real(ccs_real) :: disturbance

    dx = domain_size / real(cps)

    call get_total_num_cells(total_num_cells)
    do icell = 1, total_num_cells
      icell_global = mesh%topo%global_indices(icell)

      if (mesh%geo%x_p(1, icell) >= 0.1 * domain_size .and. mesh%geo%x_p(1, icell) <= 0.9 * domain_size .and. &
          mesh%geo%x_p(2, icell) >= 0.1 * domain_size .and. mesh%geo%x_p(2, icell) <= 0.9 * domain_size) then
        idim = 1
        disturbance = real(modulo(3 * icell_global, 17), ccs_real) / 17.0_ccs_real
        mesh%geo%x_p(idim, icell) = mesh%geo%x_p(idim, icell) + (disturbance - 0.5_ccs_real) * dx / 20.0

        idim = 2
        disturbance = real(modulo(3 * icell_global, 29), ccs_real) / 29.0_ccs_real
        mesh%geo%x_p(idim, icell) = mesh%geo%x_p(idim, icell) + (disturbance - 0.5_ccs_real) * dx / 20.0
      end if

    end do

    ! Update face interpolation
    call compute_face_interpolation(mesh)

  end subroutine

end module error_analysis
