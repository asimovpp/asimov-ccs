program test_equation_row_flush
  ! Verify that equation rows can be assembled and flushed through the public
  ! equation API. The checks cover matrix-plus-RHS rows, RHS-only rows, and
  ! automatic storage growth without losing entries already assembled.
  use testing_lib

  use constants, only: add_mode
  use fv_equations, only: scalar_transport_equation
  use kinds, only: ccs_int, ccs_real
  use mat, only: create_matrix_values, set_matrix_values_spec_ncols, set_matrix_values_spec_nrows
  use types, only: matrix_values, matrix_values_spec, vector_values
  use utils, only: clear_entries, set_mode
  use vec, only: create_vector_values

  implicit none

  real(ccs_real), parameter :: rtol = 1.0e-12_ccs_real
  real(ccs_real), parameter :: atol = 1.0e-12_ccs_real

  call init()

  call check_matrix_and_rhs_flush()
  call check_rhs_only_flush()
  call check_row_growth_preserves_entries()

  call fin()

contains

  subroutine check_matrix_and_rhs_flush()
    type(scalar_transport_equation) :: equation
    type(matrix_values) :: mat_values
    type(vector_values) :: rhs_values

    call make_work_values(3_ccs_int, mat_values, rhs_values)

    call equation%prepare_row(4_ccs_int, 3_ccs_int)
    call equation%add_row_entry(2_ccs_int, -1.0_ccs_real)
    call equation%add_row_entry(4_ccs_int, 3.5_ccs_real)
    call equation%add_row_entry(7_ccs_int, 0.25_ccs_real)
    call equation%set_rhs(-2.25_ccs_real)

    call equation%flush_row(mat_values, rhs_values)

    call assert_eq(mat_values%global_row_indices(1), 3_ccs_int, "matrix row")
    call assert_eq(mat_values%global_col_indices(1:3), [1_ccs_int, 3_ccs_int, 6_ccs_int], "matrix cols")
    call assert_close(mat_values%values(1:3), [-1.0_ccs_real, 3.5_ccs_real, 0.25_ccs_real], &
                      rtol, atol, "matrix coeffs")
    call assert_eq(rhs_values%global_indices(1), 3_ccs_int, "rhs row")
    call assert_close(rhs_values%values(1), -2.25_ccs_real, rtol, atol, "rhs value")

  end subroutine check_matrix_and_rhs_flush

  subroutine check_rhs_only_flush()
    type(scalar_transport_equation) :: equation
    type(matrix_values) :: mat_values
    type(vector_values) :: rhs_values

    call make_work_values(1_ccs_int, mat_values, rhs_values)

    call equation%prepare_row(6_ccs_int)
    call equation%set_rhs(5.75_ccs_real)

    call clear_entries(mat_values)
    call clear_entries(rhs_values)
    call equation%flush_rhs(rhs_values)

    call assert_bool(all(mat_values%global_row_indices == -1_ccs_int), "rhs-only leaves matrix rows clear")
    call assert_bool(all(mat_values%global_col_indices == -1_ccs_int), "rhs-only leaves matrix cols clear")
    call assert_eq(rhs_values%global_indices(1), 5_ccs_int, "rhs-only row")
    call assert_close(rhs_values%values(1), 5.75_ccs_real, rtol, atol, "rhs-only value")

  end subroutine check_rhs_only_flush

  subroutine check_row_growth_preserves_entries()
    type(scalar_transport_equation) :: equation
    type(matrix_values) :: mat_values
    type(vector_values) :: rhs_values

    call make_work_values(2_ccs_int, mat_values, rhs_values)

    call equation%prepare_row(3_ccs_int, 1_ccs_int)
    call equation%add_row_entry(2_ccs_int, -4.0_ccs_real)
    call equation%add_row_entry(3_ccs_int, 8.0_ccs_real)
    call equation%set_rhs(1.5_ccs_real)

    call equation%flush_row(mat_values, rhs_values)

    call assert_eq(mat_values%global_row_indices(1), 2_ccs_int, "grown row")
    call assert_eq(mat_values%global_col_indices(1:2), [1_ccs_int, 2_ccs_int], "grown cols")
    call assert_close(mat_values%values(1:2), [-4.0_ccs_real, 8.0_ccs_real], rtol, atol, "grown coeffs")
    call assert_close(rhs_values%values(1), 1.5_ccs_real, rtol, atol, "grown rhs")

  end subroutine check_row_growth_preserves_entries

  subroutine make_work_values(ncols, mat_values, rhs_values)
    integer(ccs_int), intent(in) :: ncols
    type(matrix_values), intent(out) :: mat_values
    type(vector_values), intent(out) :: rhs_values

    type(matrix_values_spec) :: mat_spec

    call set_matrix_values_spec_nrows(1_ccs_int, mat_spec)
    call set_matrix_values_spec_ncols(ncols, mat_spec)
    call create_matrix_values(mat_spec, mat_values)
    call create_vector_values(1_ccs_int, rhs_values)
    call set_mode(add_mode, mat_values)
    call set_mode(add_mode, rhs_values)

  end subroutine make_work_values

end program test_equation_row_flush
