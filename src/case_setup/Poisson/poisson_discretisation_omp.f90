
submodule (poisson_discretisation) poisson_discretisation_omp

  implicit none
  
contains

  module subroutine discretise_poisson(mesh, M)
    
    use custom_matrix, only: csr_matrix, create_new_matrix, insert_values, print_matrix
    use omp_lib

    type(ccs_mesh), intent(in) :: mesh
    class(ccs_matrix), intent(inout) :: M

    type(matrix_values_spec) :: mat_val_spec
    type(matrix_values) :: mat_coeffs
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: i, ii, j, k

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

    type(csr_matrix) :: csrm
    integer(ccs_int) :: values_per_row
    integer(ccs_int) :: value_counter
    integer(ccs_int) :: index_p
    real(ccs_real), dimension(:), allocatable :: vals
    integer(ccs_int), dimension(:), allocatable :: cols
    integer(ccs_int) :: cells_per_thread

    integer(ccs_int) :: step_by
    integer(ccs_int) :: num_threads

    !$omp parallel
      !$omp master
        num_threads = omp_get_num_threads()
        print *, "running on ", num_threads, " threads"
      !$omp end master
    !$omp end parallel

    values_per_row = 6
    cells_per_thread = 1
    step_by = cells_per_thread * num_threads
    allocate(vals(values_per_row))
    allocate(cols(values_per_row))

    ! Loop over cells
    call get_local_num_cells(mesh, local_num_cells)

    !v @todo OLD COMMENT
    !        Doing this in a loop is awful code - malloc maximum coefficients per row once,
    !        filling from front, and pass the number of coefficients to be set, requires
    !        modifying the matrix_values type and the implementation of set_values applied to
    !        matrices.
    print *, "values: ", values_per_row, step_by, local_num_cells 
    do ii = 1, local_num_cells, step_by
      print *, ii
      
      ! Will probably need a special case to handle cases where local_num_cells does not divide evenly by step_by
      call create_new_matrix(step_by, values_per_row, csrm) 
    
      !$omp parallel do &
      !$omp PRIVATE (i, &
      !$omp          loc_p, global_index_p, nnb, &
      !$omp          mat_val_spec, mat_coeffs, &
      !$omp          row, col, j, &
      !$omp          coeff_p, coeff_f, coeff_nb, &
      !$omp          loc_nb, loc_f, is_boundary, &
      !$omp          A, &
      !$omp          vals, cols, index_p, value_counter)
      do i = 1, step_by
        vals(:) = 0.0 
        cols(:) = 0 
        value_counter = 1

        index_p = ii - 1 + i
        !print *, "indices", ii, i, index_p

        call set_cell_location(mesh, index_p, loc_p)
        call get_global_index(loc_p, global_index_p)
        call count_neighbours(loc_p, nnb)

        !call set_matrix_values_spec_nrows(1_ccs_int, mat_val_spec)
        !call set_matrix_values_spec_ncols((1_ccs_int + nnb), mat_val_spec)
        !call create_matrix_values(mat_val_spec, mat_coeffs)
        !call set_mode(insert_mode, mat_coeffs)

        row = global_index_p
        coeff_p = 0.0_ccs_real
        !print *, "row, nnb", row, nnb
        row = i !XXX: temporary bodge

        ! Loop over faces
        do j = 1, nnb
          call set_neighbour_location(loc_p, j, loc_nb)
          call get_boundary_status(loc_nb, is_boundary)

          if (.not. is_boundary) then
            ! Interior face

            call set_face_location(mesh, index_p, j, loc_f)
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

          !call set_row(row, mat_coeffs)
          !call set_col(col, mat_coeffs)
          !call set_entry(coeff_nb, mat_coeffs)
          cols(value_counter) = col
          vals(value_counter) = coeff_nb
          value_counter = value_counter + 1


        end do

        ! Add the diagonal entry
        col = row
        !call set_row(row, mat_coeffs)
        !call set_col(col, mat_coeffs)
        !call set_entry(coeff_p, mat_coeffs)

        cols(value_counter) = col
        vals(value_counter) = coeff_p
        value_counter = value_counter + 1
         
        !print *, "inserted this many values", value_counter
        call insert_values(i, vals, cols, csrm)

      end do 
      !$omp end parallel do

      call print_matrix(csrm)

      ! Set the values
      ! XXX: I'm amazed this doesn't seem to cause conflicts...
      call set_matrix_values_spec_nrows(step_by, mat_val_spec)
      call set_matrix_values_spec_ncols(values_per_row, mat_val_spec)
      call create_matrix_values(mat_val_spec, mat_coeffs)
      call set_mode(insert_mode, mat_coeffs)

      do j = 1, step_by
        call set_row(ii - 1 + j, mat_coeffs)
        do k = 1, values_per_row 
          call set_col(csrm%columns(k), mat_coeffs)
          call set_entry(csrm%values(k), mat_coeffs)
        end do
      end do


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
  
end submodule poisson_discretisation_omp
