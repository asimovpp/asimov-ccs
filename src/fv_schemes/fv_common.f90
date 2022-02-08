!> @brief Submodule file fv_CDS.smod
!> @build CDS
!
!> @details An implementation of the finite volume method using the CDS scheme

submodule (fv) fv_common

  use types, only : face_locator, set_face_location
  use mesh_utils, only : face_area

  implicit none

contains

  !> @brief Computes fluxes and assign to matrix and RHS
  !
  !> @param[in] u, v      - arrays containing velocity fields in x, y directions
  !> @param[in] cell_mesh - the mesh being used
  !> @param[in] cps       - the number of cells per side in the (square) mesh
  !> @param[in,out] mat   - Data structure containing matrix to be filled
  !> @param[in,out] vec   - Data structure containing RHS vector to be filled
  module subroutine compute_fluxes(u, v, cell_mesh, cps, M, vec)
    class(field), intent(in) :: u, v
    type(mesh), intent(in) :: cell_mesh
    integer(accs_int), intent(in) :: cps
    class(matrix), intent(inout), allocatable :: M   
    class(vector), intent(inout) :: vec   

    integer(accs_int) :: n_int_cells
    
    ! Loop over cells computing advection and diffusion fluxes
    n_int_cells = calc_matrix_nnz()
    call compute_interior_coeffs(u, v, cell_mesh, n_int_cells, cps, M)

    ! Loop over boundaries
    call compute_boundary_coeffs(u, v, cell_mesh, cps, M, vec)

  end subroutine compute_fluxes

  !> @brief Returns the number of entries per row that are non-zero
  !
  !> @details Note: this assumes a square 2d grid
  !
  !> @param[out] nnz - number of non-zero entries per row
  pure function calc_matrix_nnz() result(nnz)
    integer(accs_int) :: nnz

    nnz = 5_accs_int
  end function calc_matrix_nnz

  !> @brief Computes the matrix coefficient for cells in the interior of the mesh
  !
  !> @param[in] u, v        - Field structures containing x, y velocities
  !> @param[in] cell_mesh   - Mesh structure
  !> @param[in] n_int_cells - number of cells in the interior of the mesh
  !> @param[in] cps         - number of cells per side
  !> @param[in,out] mat     - Matrix structure being assigned
  subroutine compute_interior_coeffs(u, v, cell_mesh, n_int_cells, cps, M)
    use constants, only : add_mode
    use types, only: matrix_values, cell_locator, face_locator, neighbour_locator, &
                     set_cell_location, set_face_location, set_neighbour_location
    use utils, only: pack_entries, set_values
    use mesh_utils, only:  global_index, count_neighbours, face_area, boundary_status

    class(field), intent(in) :: u, v
    type(mesh), intent(in) :: cell_mesh
    integer(accs_int), intent(in) :: n_int_cells
    integer(accs_int), intent(in) :: cps
    class(matrix), allocatable :: M

    type(matrix_values) :: mat_coeffs
    type(cell_locator) :: self_loc
    type(neighbour_locator) :: ngb_loc
    type(face_locator) :: face_loc
    integer(accs_int) :: self_idx, ngb_idx, local_idx
    integer(accs_int) :: j
    integer(accs_int) :: mat_counter
    integer(accs_int) :: n_ngb
    real(accs_real) :: face_surface_area
    real(accs_real) :: diff_coeff, diff_coeff_total
    real(accs_real) :: adv_coeff, adv_coeff_total
    logical :: is_boundary

    mat_coeffs%mode = add_mode

    allocate(mat_coeffs%rglob(1))
    allocate(mat_coeffs%cglob(n_int_cells))
    allocate(mat_coeffs%val(n_int_cells))

    do local_idx = 1, cell_mesh%nlocal
      ! Calculate contribution from neighbours
      call set_cell_location(self_loc, cell_mesh, local_idx)
      call global_index(self_loc, self_idx)
      call count_neighbours(self_loc, n_ngb)
      mat_counter = 1
      adv_coeff_total = 0.0_accs_real
      diff_coeff_total = 0.0_accs_real
      do j = 1, n_ngb
        call set_neighbour_location(ngb_loc, self_loc, j)
        call global_index(ngb_loc, ngb_idx)
        call boundary_status(ngb_loc, is_boundary)

        diff_coeff = calc_diffusion_coeff(local_idx, j, cell_mesh)
        if (.not. is_boundary) then
          call set_face_location(face_loc, cell_mesh, local_idx, j)
          call face_area(face_loc, face_surface_area)

          select type(u)
          type is(central_field)
            select type(v)
            type is(central_field)
              call calc_advection_coeff(ngb_idx, self_idx, face_surface_area, cps, u, v, 0, adv_coeff)
            class default
              print *, 'invalid velocity field discretisation'
              stop
            end select
          type is(upwind_field)
            select type(v)
            type is(upwind_field)
              call calc_advection_coeff(ngb_idx, self_idx, face_surface_area, cps, u, v, 0, adv_coeff)
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
  !> @param[in] u, v        - Field structures containing x, y velocities
  !> @param[in] cell_mesh   - Mesh structure
  !> @param[in] cps         - number of cells per side
  !> @param[in,out] M       - Matrix structure being assigned
  !> @param[in,out] b       - vector structure being assigned
  subroutine compute_boundary_coeffs(u, v, cell_mesh, cps, M, b)
    use constants, only : insert_mode, add_mode, ndim
    use types, only: matrix_values, vector_values, cell_locator, face_locator, neighbour_locator, &
                     set_cell_location, set_face_location, set_neighbour_location
    use utils, only: pack_entries, set_values
    use mesh_utils, only:  global_index, count_neighbours, face_area, &
                           boundary_status

    class(field), intent(in) :: u, v
    type(mesh), intent(in) :: cell_mesh
    integer(accs_int), intent(in) :: cps
    class(matrix), intent(inout) :: M
    class(vector), intent(inout) :: b

    type(matrix_values) :: mat_coeffs
    type(vector_values) :: b_coeffs
    type(cell_locator) :: self_loc
    type(neighbour_locator) :: ngb_loc
    type(face_locator) :: face_loc
    integer(accs_int) :: self_idx, ngb_idx, local_idx
    integer(accs_int) :: j
    integer(accs_int) :: bc_counter
    integer(accs_int) :: row, col
    integer(accs_int) :: n_ngb, mesh_ngb_idx
    real(accs_real) :: face_surface_area
    real(accs_real) :: diff_coeff
    real(accs_real) :: adv_coeff
    real(accs_real) :: BC_value
    real(accs_real) :: n_value, w_value
    logical :: is_boundary

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
    do local_idx = 1, cell_mesh%nlocal
      call set_cell_location(self_loc, cell_mesh, local_idx)
      call global_index(self_loc, self_idx)
      call count_neighbours(self_loc, n_ngb)
      ! Calculate contribution from neighbours
      do j = 1, n_ngb
        call set_neighbour_location(ngb_loc, self_loc, j)
        call global_index(ngb_loc, ngb_idx)
        call boundary_status(ngb_loc, is_boundary)

        call set_face_location(face_loc, cell_mesh, local_idx, j)
        call face_area(face_loc, face_surface_area)

        mesh_ngb_idx = cell_mesh%nbidx(j, local_idx)
        diff_coeff = calc_diffusion_coeff(local_idx, j, cell_mesh)
        if (is_boundary .and. mesh_ngb_idx == -1) then
        !if (is_boundary) then
          ! Set Dirichlet BCs at left boundary
          select type(u)
          type is(central_field)
            select type(v)
            type is(central_field)
              call calc_advection_coeff(ngb_idx, self_idx, face_surface_area, cps, u, v, mesh_ngb_idx, adv_coeff)
            class default
              print *, 'invalid velocity field discretisation'
              stop
            end select
          type is(upwind_field)
            select type(v)
            type is(upwind_field)
              call calc_advection_coeff(ngb_idx, self_idx, face_surface_area, cps, u, v, mesh_ngb_idx, adv_coeff)
            class default
              print *, 'invalid velocity field discretisation'
              stop
            end select
          class default
            print *, 'invalid velocity field discretisation'
            stop
          end select

          call calc_cell_coords(self_idx, cps, row, col)
          BC_value = -(1.0_accs_real - real(row, accs_real)/real(cps, accs_real)) * w_value
          !BC_value = compute_boundary_values(row, col, mesh_ngb_idx)
          call pack_entries(b_coeffs, 1, self_idx, (adv_coeff + diff_coeff)*BC_value)
          call pack_entries(mat_coeffs, 1, 1, self_idx, self_idx, -(adv_coeff + diff_coeff))
          call set_values(b_coeffs, b)
          call set_values(mat_coeffs, M)
          bc_counter = bc_counter + 1
        else if (is_boundary .and. mesh_ngb_idx == -4) then
          ! Set Dirichlet BCs at top boundary
          select type(u)
          type is(central_field)
            select type(v)
            type is(central_field)
              call calc_advection_coeff(ngb_idx, self_idx, face_surface_area, cps, u, v, mesh_ngb_idx, adv_coeff)
            class default
              print *, 'invalid velocity field discretisation'
              stop
            end select
          type is(upwind_field)
            select type(v)
            type is(upwind_field)
              call calc_advection_coeff(ngb_idx, self_idx, face_surface_area, cps, u, v, mesh_ngb_idx, adv_coeff)
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

  !function compute_boundary_values(row, col, BC_type) return(BC_value)
  !  integer, intent(in) :: row
  !  integer, intent(in) :: col
  !  integer, intent(in) :: BC_type
  !  real(accs_real) :: BC_value

  !  BC_value = 0.0_accs_real
  !end function compute_boundary_values

  !> @brief Sets the diffusion coefficient
  !
  !> @param[in] self_idx - the current cell index
  !> @param[in] ngb_idx  - the neigbouring cell index
  !> @param[out] coeff   - the diffusion coefficient
  module function calc_diffusion_coeff(local_self_idx, local_ngb_idx, cell_mesh) result(coeff)
    integer(accs_int), intent(in) :: local_self_idx
    integer(accs_int), intent(in) :: local_ngb_idx
    type(mesh), intent(in) :: cell_mesh
    real(accs_real) :: coeff

    type(face_locator) :: face_location
    real(accs_real) :: face_surface_area
    real(accs_real), parameter :: diffusion_factor = 1.e-2 

    call set_face_location(face_location, cell_mesh, local_self_idx, local_ngb_idx)
    call face_area(face_location, face_surface_area)

    coeff = -face_surface_area*diffusion_factor/(0.5_accs_real * cell_mesh%h)
  end function calc_diffusion_coeff

  !> @brief Calculates mass flux across given face. Note: assumes rho = 1 and uniform grid
  !
  !> @param[in] edge_len           - Edge length
  !> @param[in] u, v               - Field structures containing x, y velocities
  !> @param[in] ngb_row, ngb_col   - Row and column of given neighbour in mesh
  !> @param[in] self_row, self_col - Row and column of given cell in mesh
  !> @param[in] BC_flag            - Flag to indicate if neighbour is a boundary cell
  !> @param[out] flux              - The flux across the boundary
  module function calc_mass_flux(u, v, ngb_row, ngb_col, self_row, self_col, face_area, BC_flag) result(flux)
    class(field), intent(in) :: u, v
    integer(accs_int), intent(in) :: ngb_row, ngb_col
    integer(accs_int), intent(in) :: self_row, self_col
    real(accs_real), intent(in) :: face_area
    integer(accs_int), intent(in) :: BC_flag

    real(accs_real) :: flux

    flux = 0.0_accs_real

    if (BC_flag == 0) then
      if (ngb_col - self_col == 1) then
        flux = 0.5_accs_real*(u%val(ngb_col, ngb_row) + u%val(self_col, self_row)) * face_area
      else if (ngb_col - self_col == -1) then
        flux = -0.5_accs_real*(u%val(ngb_col, ngb_row) + u%val(self_col, self_row)) * face_area
      else if (ngb_row - self_row == 1) then
        flux = 0.5_accs_real*(v%val(ngb_col, ngb_row) + v%val(self_col, self_row)) * face_area
      else 
        flux = -0.5_accs_real*(v%val(ngb_col, ngb_row) + v%val(self_col, self_row)) * face_area
      end if
    else if (BC_flag == -1 .or. BC_flag == -2) then
      flux = u%val(self_col, self_row) * face_area
    else if (BC_flag == -3 .or. BC_flag == -4) then
      flux = v%val(self_col, self_row) * face_area
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
