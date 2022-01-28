!> @brief Submodule file fv_CDS.smod
!> @build CDS
!
!> @details An implementation of the finite volume method using the CDS scheme

submodule (fv) fv_common

  implicit none

contains

  !> @brief Computes fluxes and assign to matrix and RHS
  !
  !> @param[in,out] mat - Data structure containing matrix to be filled
  !> @param[in,out] vec - Data structure containing RHS vector to be filled
  !> @param[in] u, v - arrays containing velocity fields in x, y directions
  !> @param[in] cell_mesh - the mesh being used
  module subroutine compute_fluxes(M, vec, u, v, cell_mesh)
    class(matrix), intent(inout), allocatable :: M   
    class(vector), intent(inout) :: vec   
    class(field), intent(in) :: u, v
    type(mesh), intent(in) :: cell_mesh
    
    integer(accs_int) :: cps
    integer(accs_int) :: n_int_cells
    
    cps = int(sqrt(real(cell_mesh%n)))

    ! Loop over cells computing advection and diffusion fluxes
    n_int_cells = calc_matrix_nnz()
    call compute_interior_coeffs(M, u, v, cell_mesh, n_int_cells, cps)

    ! Loop over boundaries
    call compute_boundary_coeffs(M, vec, u, v, cell_mesh, cps)

  end subroutine compute_fluxes

  !> @brief Returns the number of entries per row that are non-zero
  !
  !> @details Note: this assumes a square 2d grid
  !
  !> @param[out] nnz - number of non-zero entries per row
  pure function calc_matrix_nnz() result(nnz)
    integer(accs_int) :: nnz

    nnz = 5
  end function calc_matrix_nnz

  !> @brief Computes the matrix coefficient for cells in the interior of the mesh
  !
  !> @param[in,out] mat     - Matrix structure being assigned
  !> @param[in] u, v        - Field structures containing x, y velocities
  !> @param[in] cell_mesh   - Mesh structure
  !> @param[in] n_int_cells - number of cells in the interior of the mesh
  !> @param[in] cps         - number of cells per side
  subroutine compute_interior_coeffs(M, u, v, cell_mesh, n_int_cells, cps)
    use constants, only : add_mode
    use types, only: matrix_values
    use utils, only: pack_entries, set_values

    class(matrix), allocatable :: M
    class(field), intent(in) :: u, v
    type(mesh), intent(in) :: cell_mesh
    integer(accs_int), intent(in) :: n_int_cells
    integer(accs_int), intent(in) :: cps

    type(matrix_values) :: mat_coeffs
    integer(accs_int) :: self_idx, ngb_idx, local_idx
    integer(accs_int) :: j
    real(accs_real) :: face_area
    real(accs_real) :: diff_coeff, diff_coeff_total
    real(accs_real) :: adv_coeff, adv_coeff_total
    integer(accs_int) :: mat_counter

    mat_coeffs%mode = add_mode

    allocate(mat_coeffs%rglob(1))
    allocate(mat_coeffs%cglob(n_int_cells))
    allocate(mat_coeffs%val(n_int_cells))

    do local_idx = 1, cell_mesh%nlocal
      ! Calculate contribution from neighbours
      self_idx = cell_mesh%idx_global(local_idx)
      mat_counter = 1
      adv_coeff_total = 0.0_accs_real
      diff_coeff_total = 0.0_accs_real
      do j = 1, cell_mesh%nnb(local_idx)
        ngb_idx = cell_mesh%nbidx(j, local_idx)
        face_area = cell_mesh%Af(j, local_idx)
        diff_coeff = calc_diffusion_coeff()
        if (ngb_idx > 0) then
          select type(u)
          type is(central_field)
            select type(v)
            type is(central_field)
              call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, cps, u, v, 0)
            class default
              print *, 'invalid velocity field discretisation'
              stop
            end select
          type is(upwind_field)
            select type(v)
            type is(upwind_field)
              call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, cps, u, v, 0)
            class default
              print *, 'invalid velocity field discretisation'
              stop
            end select
          class default
            print *, 'invalid velocity field discretisation'
            stop
          end select

          call pack_entries(mat_coeffs, 1, mat_counter, self_idx, ngb_idx, adv_coeff + diff_coeff)
          mat_counter = mat_counter + 1
          adv_coeff_total = adv_coeff_total + adv_coeff
          diff_coeff_total = diff_coeff_total + diff_coeff
        else
          call pack_entries(mat_coeffs, 1, mat_counter, self_idx, -1, 0.0_accs_real)
          mat_counter = mat_counter + 1
        end if
      end do
      call pack_entries(mat_coeffs, 1, mat_counter, self_idx, self_idx, -(adv_coeff_total + diff_coeff_total))
      mat_counter = mat_counter + 1
      call set_values(mat_coeffs, M)
    end do

    deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)
  end subroutine compute_interior_coeffs

  !> @brief Computes the matrix coefficient for cells on the boundary of the mesh
  !
  !> @param[in,out] mat     - Matrix structure being assigned
  !> @param[in] u, v        - Field structures containing x, y velocities
  !> @param[in] cell_mesh   - Mesh structure
  !> @param[in] cps         - number of cells per side
  subroutine compute_boundary_coeffs(M, b, u, v, cell_mesh, cps)
    use constants, only : insert_mode, add_mode
    use types, only: matrix_values, vector_values
    use utils, only: pack_entries, set_values

    class(matrix), intent(inout) :: M
    class(vector), intent(inout) :: b
    class(field), intent(in) :: u, v
    type(mesh), intent(in) :: cell_mesh
    integer(accs_int), intent(in) :: cps

    type(matrix_values) :: mat_coeffs
    type(vector_values) :: b_coeffs

    integer(accs_int) :: self_idx, ngb_idx, local_idx
    integer(accs_int) :: j
    integer(accs_int) :: bc_counter
    integer(accs_int) :: row, col
    real(accs_real) :: face_area
    real(accs_real) :: diff_coeff
    real(accs_real) :: adv_coeff
    real(accs_real) :: BC_value
    real(accs_real) :: n_value, w_value

    mat_coeffs%mode = add_mode
    b_coeffs%mode = add_mode

    allocate(mat_coeffs%rglob(1))
    allocate(mat_coeffs%cglob(1))
    allocate(mat_coeffs%val(1))
    allocate(b_coeffs%idx(1))
    allocate(b_coeffs%val(1))

    n_value = 0.0_accs_real
    w_value = 1.0_accs_real
    bc_counter = 1
    diff_coeff = calc_diffusion_coeff()
    do local_idx = 1, cell_mesh%nlocal
      self_idx = cell_mesh%idx_global(local_idx)
      ! Calculate contribution from neighbours
      do j = 1, cell_mesh%nnb(local_idx)
        ngb_idx = cell_mesh%nbidx(j, local_idx)
        face_area = cell_mesh%Af(j, local_idx)
        if (ngb_idx == -1) then
          ! Set Dirichlet BCs at left boundary
          select type(u)
          type is(central_field)
            select type(v)
            type is(central_field)
              call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, cps, u, v, ngb_idx)
            class default
              print *, 'invalid velocity field discretisation'
              stop
            end select
          type is(upwind_field)
            select type(v)
            type is(upwind_field)
              call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, cps, u, v, ngb_idx)
            class default
              print *, 'invalid velocity field discretisation'
              stop
            end select
          class default
            print *, 'invalid velocity field discretisation'
            stop
          end select

          call calc_cell_coords(self_idx, cps, row, col)
          BC_value = -(1.0_accs_real - real(row, kind=accs_real)/real(cps, kind=accs_real)) * w_value
          call pack_entries(b_coeffs, 1, self_idx, (adv_coeff + diff_coeff)*BC_value)
          call pack_entries(mat_coeffs, 1, 1, self_idx, self_idx, -(adv_coeff + diff_coeff))
          call set_values(b_coeffs, b)
          call set_values(mat_coeffs, M)
          bc_counter = bc_counter + 1
        else if (ngb_idx == -4) then
          ! Set Dirichlet BCs at top boundary
          select type(u)
          type is(central_field)
            select type(v)
            type is(central_field)
              call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, cps, u, v, ngb_idx)
            class default
              print *, 'invalid velocity field discretisation'
              stop
            end select
          type is(upwind_field)
            select type(v)
            type is(upwind_field)
              call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, cps, u, v, ngb_idx)
            class default
              print *, 'invalid velocity field discretisation'
              stop
            end select
          class default
            print *, 'invalid velocity field discretisation'
            stop
          end select

          call pack_entries(b_coeffs, 1, self_idx, (adv_coeff + diff_coeff)*n_value)
          call pack_entries(mat_coeffs, 1, 1, self_idx, self_idx, -(adv_coeff + diff_coeff))
          call set_values(b_coeffs, b)
          call set_values(mat_coeffs, M)
          bc_counter = bc_counter + 1
        end if
      end do
    end do
    deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)
    deallocate(b_coeffs%idx, b_coeffs%val)
  end subroutine compute_boundary_coeffs

  !> @brief Sets the diffusion coefficient
  !
  !> @param[out] coeff - the diffusion coefficient
  pure function calc_diffusion_coeff() result(coeff)
    real(accs_real) :: coeff

    real(accs_real), parameter :: diffusion_factor = 1.e-2

    coeff = -2.*diffusion_factor
  end function calc_diffusion_coeff

  !> @brief Calculates mass flux across given face. Note: assumes rho = 1 and uniform grid
  !
  !> @param[in] edge_len           - Edge length
  !> @param[in] u, v               - Field structures containing x, y velocities
  !> @param[in] ngb_row, ngb_col   - Row and column of given neighbour in mesh
  !> @param[in] self_row, self_col - Row and column of given cell in mesh
  !> @param[in] BC_flag            - Flag to indicate if neighbour is a boundary cell
  !> @param[out] flux              - The flux across the boundary
  module function calc_mass_flux(edge_len, u, v, ngb_row, ngb_col, self_row, self_col, BC_flag) result(flux)
    real(accs_real), intent(in) :: edge_len
    class(field), intent(in) :: u, v
    integer(accs_int), intent(in) :: ngb_row, ngb_col
    integer(accs_int), intent(in) :: self_row, self_col
    integer(accs_int), intent(in) :: BC_flag

    real(accs_real) :: flux

    flux = 0.

    if (BC_flag == 0) then
      if (ngb_col - self_col == 1) then
        flux = 0.25*(u%val(ngb_col, ngb_row) + u%val(self_col, self_row)) * edge_len
      else if (ngb_col - self_col == -1) then
        flux = -0.25*(u%val(ngb_col, ngb_row) + u%val(self_col, self_row)) * edge_len
      else if (ngb_row - self_row == 1) then
        flux = 0.25*(v%val(ngb_col, ngb_row) + v%val(self_col, self_row)) * edge_len
      else 
        flux = -0.25*(v%val(ngb_col, ngb_row) + v%val(self_col, self_row)) * edge_len
      end if
    else if (BC_flag == -1 .or. BC_flag == -2) then
      flux = u%val(self_col, self_row) * edge_len
    else if (BC_flag == -3 .or. BC_flag == -4) then
      flux = v%val(self_col, self_row) * edge_len
    else
      print *, 'invalid BC flag'
      stop
    end if
  end function calc_mass_flux

  !> @brief Calculates the row and column indices from flattened vector index. Assumes square mesh
  !
  !> @param[in] idx  - cell index
  !> @param[in] cps  - number of cells per side
  !> @param[out] row - cell row within mesh
  !> @param[out] col - cell column within mesh
  module subroutine calc_cell_coords(idx, cps, row, col)
    integer(accs_int), intent(in) :: idx, cps
    integer(accs_int), intent(out) :: row, col

    col = modulo(idx-1,cps) + 1 
    row = (idx-1)/cps + 1
  end subroutine calc_cell_coords
  

end submodule fv_common
