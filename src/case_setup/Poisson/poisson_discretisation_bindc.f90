
submodule (poisson_discretisation) poisson_discretisation_bindc

  implicit none
  
contains

  module subroutine discretise_poisson(mesh, M)

    type(ccs_mesh), intent(in) :: mesh
    class(ccs_matrix), intent(inout) :: M

    type(matrix_values_spec) :: mat_val_spec
    type(matrix_values) :: mat_coeffs
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: i, j

    integer(ccs_int) :: row, col
    real(ccs_real) :: coeff_f, coeff_p, coeff_nb

    type(face_locator) :: loc_f
    real(ccs_real) :: A

    integer(ccs_int) :: global_index_p
    type(cell_locator) :: loc_p
    integer(ccs_int) :: nnb

    type(neighbour_locator) :: loc_nb
    logical :: is_boundary
    integer(ccs_int) :: global_index_nb

    ! Loop over cells
    call get_local_num_cells(mesh, local_num_cells)
    do i = 1, local_num_cells
      !^ @todo Doing this in a loop is awful code - malloc maximum coefficients per row once,
      !        filling from front, and pass the number of coefficients to be set, requires
      !        modifying the matrix_values type and the implementation of set_values applied to
      !        matrices.
      call set_cell_location(mesh, i, loc_p)
      call get_global_index(loc_p, global_index_p)
      call count_neighbours(loc_p, nnb)

      call set_matrix_values_spec_nrows(1_ccs_int, mat_val_spec)
      call set_matrix_values_spec_ncols((1_ccs_int + nnb), mat_val_spec)
      call create_matrix_values(mat_val_spec, mat_coeffs)
      call set_mode(insert_mode, mat_coeffs)

      row = global_index_p
      coeff_p = 0.0_ccs_real

      ! Loop over faces
      do j = 1, nnb
        call set_neighbour_location(loc_p, j, loc_nb)
        call get_boundary_status(loc_nb, is_boundary)

        if (.not. is_boundary) then
          ! Interior face

          call set_face_location(mesh, i, j, loc_f)
          call get_face_area(loc_f, A)
          coeff_f = (1.0 / mesh%geo%h) * A

          call get_global_index(loc_nb, global_index_nb)

          coeff_p = coeff_p - coeff_f
          coeff_nb = coeff_f
          col = global_index_nb
        else
          col = -1
          coeff_nb = 0.0_ccs_real
        end if

        call set_row(row, mat_coeffs)
        call set_col(col, mat_coeffs)
        call set_entry(coeff_nb, mat_coeffs)

      end do

      ! Add the diagonal entry
      col = row
      call set_row(row, mat_coeffs)
      call set_col(col, mat_coeffs)
      call set_entry(coeff_p, mat_coeffs)

      ! Set the values
      call set_values(mat_coeffs, M)

    end do

  end subroutine discretise_poisson

  module subroutine apply_dirichlet_bcs(mesh, M, b)

    implicit none

    type(ccs_mesh), intent(in) :: mesh
    class(ccs_matrix), intent(inout) :: M
    class(ccs_vector), intent(inout) :: b

    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: i, j
    real(ccs_real) :: boundary_coeff, boundary_val

    integer(ccs_int) :: idx, row, col
    real(ccs_real) :: r, coeff

    type(vector_values) :: vec_values
    type(matrix_values_spec) :: mat_val_spec
    type(matrix_values) :: mat_coeffs

    type(face_locator) :: loc_f
    real(ccs_real) :: A

    type(cell_locator) :: loc_p
    integer(ccs_int) :: global_index_p
    type(neighbour_locator) :: loc_nb

    integer(ccs_int) :: nnb
    logical :: is_boundary

    integer(ccs_int) :: nrows_working_set

    call set_matrix_values_spec_nrows(1_ccs_int, mat_val_spec)
    call set_matrix_values_spec_ncols(1_ccs_int, mat_val_spec)
    call create_matrix_values(mat_val_spec, mat_coeffs)
    call set_mode(add_mode, mat_coeffs)

    nrows_working_set = 1_ccs_int
    call create_vector_values(nrows_working_set, vec_values)
    call set_mode(add_mode, vec_values)

    call get_local_num_cells(mesh, local_num_cells)
    do i = 1, local_num_cells
      if (minval(mesh%topo%nb_indices(:, i)) < 0) then
        call clear_entries(mat_coeffs)
        call clear_entries(vec_values)
        call set_cell_location(mesh, i, loc_p)
        call get_global_index(loc_p, global_index_p)
        coeff = 0.0_ccs_real
        r = 0.0_ccs_real

        row = global_index_p
        col = global_index_p
        idx = global_index_p

        call count_neighbours(loc_p, nnb)
        do j = 1, nnb

          call set_neighbour_location(loc_p, j, loc_nb)
          call get_boundary_status(loc_nb, is_boundary)

          if (is_boundary) then
            call set_face_location(mesh, i, j, loc_f)
            call get_face_area(loc_f, A)
            boundary_coeff = (2.0 / mesh%geo%h) * A
            boundary_val = eval_solution(loc_f)

            ! Coefficient
            coeff = coeff - boundary_coeff

            ! RHS vector
            r = r - boundary_val * boundary_coeff
          end if

        end do

        call set_row(row, mat_coeffs)
        call set_col(col, mat_coeffs)
        call set_entry(coeff, mat_coeffs)

        call set_row(row, vec_values)
        call set_entry(r, vec_values)

        call set_values(mat_coeffs, M)
        call set_values(vec_values, b)

      end if
    end do

    deallocate (vec_values%global_indices)
    deallocate (vec_values%values)

  end subroutine apply_dirichlet_bcs
  
end submodule poisson_discretisation_bindc
