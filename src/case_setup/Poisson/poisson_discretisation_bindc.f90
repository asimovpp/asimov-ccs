
submodule (poisson_discretisation) poisson_discretisation_bindc

  use iso_c_binding

  use meshing, only: get_local_index

  implicit none

  interface
    subroutine discretise_poisson_kernel(nrows, nnz_pr, h, &
         mesh_neighbours, mesh_face_areas, csr_values) bind(c)
      use iso_c_binding

      integer(c_int) :: nrows, nnz_pr
      real(c_double) :: h
      integer(c_int), dimension(nrows * nnz_pr) :: mesh_neighbours
      real(c_double), dimension(nrows * nnz_pr) :: mesh_face_areas
      real(c_double), dimension(nrows * (nnz_pr + 1)) :: csr_values
    end subroutine discretise_poisson_kernel
  end interface

  interface
    subroutine apply_dirichlet_bcs_kernel(nrows, nnz_pr, h, &
          mesh_neighbours, mesh_face_areas, mesh_face_yloc, diag_values, rhs_values) bind(c)
      use iso_c_binding

      integer(c_int) :: nrows, nnz_pr
      real(c_double) :: h
      integer(c_int), dimension(nrows * nnz_pr) :: mesh_neighbours
      real(c_double), dimension(nrows * nnz_pr) :: mesh_face_areas
      real(c_double), dimension(nrows * nnz_pr) :: mesh_face_yloc
      real(c_double), dimension(nrows) :: diag_values
      real(c_double), dimension(nrows) :: rhs_values
    end subroutine apply_dirichlet_bcs_kernel
  end interface

contains

  module subroutine discretise_poisson(mesh, M)

    type(ccs_mesh), intent(in) :: mesh
    class(ccs_matrix), intent(inout) :: M

    type(matrix_values_spec) :: mat_val_spec
    type(matrix_values) :: mat_coeffs
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: i, j

    integer(ccs_int) :: row, col
    real(ccs_real) :: coeff_p, coeff_nb

    integer(ccs_int) :: global_index_p
    type(cell_locator) :: loc_p
    integer(ccs_int) :: nnb

    integer(ccs_int) :: global_index_nb

    integer(ccs_int), dimension(:), allocatable :: flat_neighbours
    real(ccs_real), dimension(:), allocatable :: flat_areas

    integer(ccs_int) :: index_nb
    real(ccs_real), dimension(:), allocatable :: csr_values
    integer(ccs_int) :: csr_idx, face_idx

    ! Fake a flat mesh data structure
    call get_local_num_cells(mesh, local_num_cells)
    call flatten_mesh(mesh, flat_neighbours, flat_areas)

    allocate(csr_values((4 + 1) * local_num_cells)) ! 1 entry per neighbour + diagonal coefficient

    call discretise_poisson_kernel(local_num_cells, 4, mesh%geo%h, &
         flat_neighbours, flat_areas, csr_values)

    ! Copy CSR-like structure into PETSc
    do i = 1, local_num_cells
      !^ @todo Doing this in a loop is awful code - malloc maximum coefficients per row once,
      !        filling from front, and pass the number of coefficients to be set, requires
      !        modifying the matrix_values type and the implementation of set_values applied to
      !        matrices.
      call set_cell_location(mesh, i, loc_p)
      call get_global_index(loc_p, global_index_p)
      call count_neighbours(loc_p, nnb)
      if (nnb /= 4) then
        print *, "Hackathon assumption of 4 face per cell is broken!"
        stop 1
      end if

      call set_matrix_values_spec_nrows(1_ccs_int, mat_val_spec)
      call set_matrix_values_spec_ncols((1_ccs_int + nnb), mat_val_spec)
      call create_matrix_values(mat_val_spec, mat_coeffs)
      call set_mode(insert_mode, mat_coeffs)

      row = global_index_p
      coeff_p = 0.0_ccs_real

      ! Loop over faces
      do j = 1, nnb
        ! call set_neighbour_location(loc_p, j, loc_nb)
        ! call get_boundary_status(loc_nb, is_boundary)

        face_idx = 4 * (i - 1) + j
        csr_idx = (4 + 1) * (i - 1) + j

        index_nb = flat_neighbours(face_idx)
        ! if (.not. is_boundary) then
        if (index_nb > 0) then
          ! Interior face
          ! call get_global_index(loc_nb, global_index_nb)
          global_index_nb = mesh%topo%global_indices(index_nb)
          col = global_index_nb

          coeff_nb = csr_values(csr_idx)
        else
          col = -1
          coeff_nb = 0.0_ccs_real
        end if

        call set_row(row, mat_coeffs)
        call set_col(col, mat_coeffs)
        call set_entry(coeff_nb, mat_coeffs)

      end do
      csr_idx = ((4 + 1) * (i - 1) + 4) + 1 ! XXX: assumes 4 neighbours
      coeff_p = csr_values(csr_idx)

      ! Add the diagonal entry
      col = row
      call set_row(row, mat_coeffs)
      call set_col(col, mat_coeffs)
      call set_entry(coeff_p, mat_coeffs)

      ! Set the values
      call set_values(mat_coeffs, M)

    end do

    deallocate(csr_values)
    deallocate(flat_neighbours)
    deallocate(flat_areas)

  end subroutine discretise_poisson

  module subroutine apply_dirichlet_bcs(mesh, M, b)

    implicit none

    type(ccs_mesh), intent(in) :: mesh
    class(ccs_matrix), intent(inout) :: M
    class(ccs_vector), intent(inout) :: b

    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: i

    integer(ccs_int) :: row, col
    real(ccs_real) :: r, coeff

    type(vector_values) :: vec_values
    type(matrix_values_spec) :: mat_val_spec
    type(matrix_values) :: mat_coeffs

    type(cell_locator) :: loc_p
    integer(ccs_int) :: global_index_p

    integer(ccs_int), dimension(:), allocatable :: flat_neighbours
    real(ccs_real), dimension(:), allocatable :: flat_areas
    real(ccs_real), dimension(:), allocatable :: flat_face_yloc
    real(ccs_real), dimension(:), allocatable :: diag_values
    real(ccs_real), dimension(:), allocatable :: rhs_values

    integer(ccs_int) :: nrows_working_set


    call get_local_num_cells(mesh, local_num_cells)

    call flatten_mesh(mesh, flat_neighbours, flat_areas, flat_face_yloc)

    allocate(diag_values(local_num_cells))
    allocate(rhs_values(local_num_cells))
    diag_values(:) = 0.0_ccs_real
    rhs_values(:) = 0.0_ccs_real

    call apply_dirichlet_bcs_kernel(local_num_cells, 4, mesh%geo%h, &
          flat_neighbours, flat_areas, flat_face_yloc, diag_values, rhs_values)

    nrows_working_set = 1_ccs_int
    call create_vector_values(nrows_working_set, vec_values)
    call set_mode(add_mode, vec_values)

 
    call set_matrix_values_spec_nrows(1_ccs_int, mat_val_spec)
    call set_matrix_values_spec_ncols(1_ccs_int, mat_val_spec)
    call create_matrix_values(mat_val_spec, mat_coeffs)
    call set_mode(add_mode, mat_coeffs)

    call get_local_num_cells(mesh, local_num_cells)
    do i = 1, local_num_cells
      if (minval(mesh%topo%nb_indices(:, i)) < 0) then
        call clear_entries(mat_coeffs)
        call clear_entries(vec_values)
        call set_cell_location(mesh, i, loc_p)
        call get_global_index(loc_p, global_index_p)

        row = global_index_p
        col = global_index_p

        coeff = diag_values(i)
        call set_row(row, mat_coeffs)
        call set_col(col, mat_coeffs)
        call set_entry(coeff, mat_coeffs)

        r = rhs_values(i)
        call set_row(row, vec_values)
        call set_entry(r, vec_values)

        call set_values(mat_coeffs, M)
        call set_values(vec_values, b)

      end if
    end do

    deallocate (vec_values%global_indices)
    deallocate (vec_values%values)

  end subroutine apply_dirichlet_bcs

  ! Fake a flat mesh data structure.
  subroutine flatten_mesh(mesh, flat_neighbours, flat_areas, flat_face_yloc)

    type(ccs_mesh), intent(in) :: mesh
    integer(ccs_int), dimension(:), allocatable, intent(out) :: flat_neighbours
    real(ccs_real), dimension(:), allocatable, intent(out) :: flat_areas
    real(ccs_real), dimension(:), allocatable, optional, intent(out) :: flat_face_yloc

    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: index_p
    integer(ccs_int) :: index_nb

    integer(ccs_int) :: j, ctr
    real(ccs_real), dimension(ndim) :: x

    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb
    type(face_locator) :: loc_f

    real(ccs_real) :: A

    call get_local_num_cells(mesh, local_num_cells)

    allocate(flat_neighbours(4 * local_num_cells))
    allocate(flat_areas(4 * local_num_cells))
    if (present(flat_face_yloc)) then
      allocate(flat_face_yloc(4 * local_num_cells))
    end if

    ctr = 1
    do index_p = 1, local_num_cells
      call set_cell_location(mesh, index_p, loc_p)

      do j = 1, 4 ! XXX: Don't leave hardcoded!
        call set_neighbour_location(loc_p, j, loc_nb)

        call get_local_index(loc_nb, index_nb)
        flat_neighbours(ctr) = index_nb

        call set_face_location(mesh, index_p, j, loc_f)
        call get_face_area(loc_f, A)
        flat_areas(ctr) = A

        if (present(flat_face_yloc)) then
          call get_centre(loc_f, x)
          flat_face_yloc(ctr) = x(2)
        end if

        ctr = ctr + 1
      end do
    end do

  end subroutine flatten_mesh

end submodule poisson_discretisation_bindc
