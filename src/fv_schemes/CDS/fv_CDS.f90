!> @brief Submodule file fv_CDS.smod
!> @build CDS
!
!> @details An implementation of the finite volume method using the CDS scheme

submodule (fv) fv_CDS

  implicit none

contains

  module subroutine compute_fluxes(mat, vec, u, v, cell_mesh)
    use constants, only : insert_mode, add_mode
    use types, only : matrix_values, vector_values
    use utils, only : set_values, pack_entries

    class(matrix), intent(inout) :: mat   
    class(vector), intent(inout) :: vec   
    real(accs_real), dimension(:,:), intent(in) :: u, v
    type(mesh), intent(in) :: cell_mesh
    
    integer(accs_int) :: cps
    integer(accs_int) :: n_bc_cells, n_int_cells
    
    cps = int(sqrt(real(cell_mesh%n)))

    ! Loop over cells computing advection and diffusion fluxes
    n_int_cells = calc_matrix_nnz()
    call compute_interior_coeffs(mat, "CDS", u, v, cell_mesh, n_int_cells, cps)

    ! Loop over boundaries
    n_bc_cells = calc_rhs_nnz(cps)
    call compute_boundary_coeffs(mat, vec, "CDS", u, v, cell_mesh, n_bc_cells, cps)

  end subroutine compute_fluxes

  ! Note: this assumes a 2d grid
  pure function calc_matrix_nnz() result(nnz)
    integer(accs_int) :: nnz

    nnz = 5
  end function calc_matrix_nnz

  ! Note: this assumes a 2d grid
  pure function calc_rhs_nnz(cps) result(nnz)
    implicit none
    integer(accs_int), intent(in) :: cps
    integer(accs_int) :: nnz

    nnz = 2*cps
  end function calc_rhs_nnz

  subroutine compute_interior_coeffs(mat, discretisation, u, v, cell_mesh, n_int_cells, cps)
    use constants, only : insert_mode, add_mode
    use types, only: matrix_values
    use utils, only: pack_entries, set_values

    class(matrix), intent(inout) :: mat
    character(len=3), intent(in) :: discretisation
    real(accs_real), dimension(:,:), intent(in) :: u, v
    type(mesh), intent(in) :: cell_mesh
    integer(accs_int), intent(in) :: n_int_cells
    integer(accs_int), intent(in) :: cps

    type(matrix_values) :: mat_coeffs
    integer(accs_int) :: self_idx, ngb_idx
    integer(accs_int) :: j
    real(accs_real) :: face_area
    real(accs_real) :: diff_coeff, diff_coeff_total
    real(accs_real) :: adv_coeff, adv_coeff_total
    integer(accs_int) :: mat_counter

    mat_coeffs%mode = add_mode

    allocate(mat_coeffs%rglob(1))
    allocate(mat_coeffs%cglob(n_int_cells))
    allocate(mat_coeffs%val(n_int_cells))

    do self_idx = 1, cell_mesh%n
      ! Calculate contribution from neighbours
      mat_counter = 1
      adv_coeff_total = 0.0_accs_real
      diff_coeff_total = 0.0_accs_real
      do j = 1, cell_mesh%nnb(self_idx)
        ngb_idx = cell_mesh%nbidx(j, self_idx)
        face_area = cell_mesh%Af(j, self_idx)
        call calc_diffusion_coeff(diff_coeff)
        if (ngb_idx > 0) then
          call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, discretisation, cps, u, v, 0)
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
      call set_values(mat_coeffs, mat)
    end do

    deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)
  end subroutine compute_interior_coeffs

  subroutine compute_boundary_coeffs(M, b, discretisation, u, v, cell_mesh, n_bc_cells, cps)
    use constants, only : insert_mode, add_mode
    use types, only: matrix_values, vector_values
    use utils, only: pack_entries, set_values

    class(matrix), intent(inout) :: M
    class(vector), intent(inout) :: b
    character(len=3), intent(in) :: discretisation
    real(accs_real), dimension(:,:), intent(in) :: u, v
    type(mesh), intent(in) :: cell_mesh
    integer(accs_int), intent(in) :: n_bc_cells
    integer(accs_int), intent(in) :: cps

    type(matrix_values) :: mat_coeffs
    type(vector_values) :: b_coeffs

    integer(accs_int) :: self_idx, ngb_idx
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

    allocate(mat_coeffs%rglob(n_bc_cells))
    allocate(mat_coeffs%cglob(1))
    allocate(mat_coeffs%val(n_bc_cells))
    allocate(b_coeffs%idx(n_bc_cells))
    allocate(b_coeffs%val(n_bc_cells))

    n_value = 0.0_accs_real
    w_value = 1.0_accs_real
    bc_counter = 1
    call calc_diffusion_coeff(diff_coeff)
    do self_idx = 1, cell_mesh%n
      ! Calculate contribution from neighbours
      do j = 1, cell_mesh%nnb(self_idx)
        ngb_idx = cell_mesh%nbidx(j, self_idx)
        face_area = cell_mesh%Af(j, self_idx)
        if (ngb_idx == -1) then
          ! Set Dirichlet BCs at left boundary
          call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, discretisation, cps, u, v, ngb_idx)
          call calc_cell_coords(self_idx, cps, row, col)
          BC_value = -(1.0_accs_real - real(row, kind=accs_real)/real(cps, kind=accs_real)) * w_value
          call pack_entries(b_coeffs, bc_counter, self_idx, (adv_coeff + diff_coeff)*BC_value)
          call pack_entries(mat_coeffs, bc_counter, 1, self_idx, self_idx, -(adv_coeff + diff_coeff))
          bc_counter = bc_counter + 1
        else if (ngb_idx == -4) then
          ! Set Dirichlet BCs at top boundary
          call calc_advection_coeff(ngb_idx, self_idx, face_area, adv_coeff, discretisation, cps, u, v, ngb_idx)
          call pack_entries(b_coeffs, bc_counter, self_idx, (adv_coeff + diff_coeff)*n_value)
          call pack_entries(mat_coeffs, bc_counter, 1, self_idx, self_idx, -(adv_coeff + diff_coeff))
          bc_counter = bc_counter + 1
        end if
      end do
    end do
    call set_values(b_coeffs, b)
    call set_values(mat_coeffs, M)
    deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)
    deallocate(b_coeffs%idx, b_coeffs%val)
  end subroutine compute_boundary_coeffs

  ! Calculates advection coefficient for neighbouring cell 
  subroutine calc_advection_coeff(ngb_idx, self_idx, face_area, coeff, discretisation, cps, u, v, BC)
    integer(accs_int), intent(in) :: ngb_idx, self_idx
    real(accs_real), intent(in) :: face_area
    real(accs_real), intent(inout) :: coeff
    character(len=3), intent(in) :: discretisation
    integer(accs_int), intent(in) :: cps
    real(accs_real), dimension(:,:) :: u, v
    integer(accs_int), intent(in) :: BC

    integer(accs_int) :: ngb_row, ngb_col       ! neighbour coordinates within grid
    integer(accs_int) :: self_row, self_col     ! cell coordinates within grid

    ! Find where we are in the grid first
    call calc_cell_coords(ngb_idx, cps, ngb_row, ngb_col)
    call calc_cell_coords(self_idx, cps, self_row, self_col)

    coeff = calc_mass_flux(face_area, u, v, ngb_row, ngb_col, self_row, self_col, BC)
    if (discretisation == "UDS") then
      coeff = min(coeff, 0.0_accs_real)
    end if
  end subroutine calc_advection_coeff

  ! Sets diffusion coefficient
  subroutine calc_diffusion_coeff(coeff)
    real(accs_real), intent(inout) :: coeff

    real(accs_real), parameter :: diffusion_factor = 1.e-2

    coeff = -2.*diffusion_factor
  end subroutine calc_diffusion_coeff

  ! Calculates mass flux across given edge. Note: assuming rho = 1 and uniform grid
  function calc_mass_flux(edge_len, u, v, ngb_row, ngb_col, self_row, self_col, BC_flag) result(flux)
    real(accs_real), intent(in) :: edge_len
    real(accs_real), dimension(:,:), intent(in) :: u, v
    integer(accs_int), intent(in) :: ngb_row, ngb_col
    integer(accs_int), intent(in) :: self_row, self_col
    integer(accs_int), intent(in) :: BC_flag

    real(accs_real) :: flux

    flux = 0.

    if (BC_flag == 0) then
      if (ngb_col - self_col == 1) then
        flux = 0.25*(u(ngb_col, ngb_row) + u(self_col, self_row)) * edge_len
      else if (ngb_col - self_col == -1) then
        flux = -0.25*(u(ngb_col, ngb_row) + u(self_col, self_row)) * edge_len
      else if (ngb_row - self_row == 1) then
        flux = 0.25*(v(ngb_col, ngb_row) + v(self_col, self_row)) * edge_len
      else 
        flux = -0.25*(v(ngb_col, ngb_row) + v(self_col, self_row)) * edge_len
      end if
    else if (BC_flag == -1 .or. BC_flag == -2) then
      flux = u(self_col, self_row) * edge_len
    else if (BC_flag == -3 .or. BC_flag == -4) then
      flux = v(self_col, self_row) * edge_len
    else
      print *, 'invalid BC flag'
      stop
    end if
  end function calc_mass_flux

  ! Assigns source vector
  ! Calculates the row and column indices from flattened vector index
  ! Note: assumes square mesh
  subroutine calc_cell_coords(idx, cps, row, col)
    integer(accs_int), intent(in) :: idx, cps
    integer(accs_int), intent(out) :: row, col

    col = modulo(idx-1,cps) + 1 
    row = (idx-1)/cps + 1
  end subroutine calc_cell_coords
  

end submodule fv_CDS