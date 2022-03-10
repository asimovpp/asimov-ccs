!> @brief Submodule file fv_CDS.smod
!> @build CDS
!
!> @details An implementation of the finite volume method using the CDS scheme

submodule (fv) fv_common

  use types, only : face_locator
  use meshing, only : set_face_location, get_face_area, get_face_normal
  use vec, only: get_vector_data, restore_vector_data

  implicit none

contains

  !> @brief Computes fluxes and assign to matrix and RHS
  !
  !> @param[in] phi       - scalar field structure
  !> @param[in] u, v      - field structures in x, y directions
  !> @param[in] cell_mesh - the mesh being used
  !> @param[in] bcs       - the boundary conditions structure being used
  !> @param[in] cps       - the number of cells per side in the (square) mesh
  !> @param[in,out] M     - Data structure containing matrix to be filled
  !> @param[in,out] vec   - Data structure containing RHS vector to be filled
  module subroutine compute_fluxes(phi, u, v, cell_mesh, bcs, cps, M, vec)
    use petsctypes, only: vector_petsc
    class(field), intent(in) :: phi
    class(field), intent(in) :: u, v
    type(mesh), intent(in) :: cell_mesh
    type(bc_config), intent(in) :: bcs
    integer(accs_int), intent(in) :: cps
    class(matrix), allocatable, intent(inout) :: M   
    class(vector), intent(inout) :: vec   

    integer(accs_int) :: n_int_cells
    real(accs_real), dimension(:), pointer :: u_data, v_data

    associate (u_vec => u%vec, v_vec => v%vec)
      call get_vector_data(u_vec, u_data)
      call get_vector_data(v_vec, v_data)

      ! Loop over cells computing advection and diffusion fluxes
      n_int_cells = calc_matrix_nnz()
      call compute_interior_coeffs(phi, u_data, v_data, cell_mesh, n_int_cells, M)

      ! Loop over boundaries
      call compute_boundary_coeffs(phi, u_data, v_data, cell_mesh, bcs, cps, M, vec)
      
      call restore_vector_data(u_vec, u_data)
      call restore_vector_data(v_vec, v_data)
    end associate

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
  !> @param[in] phi         - scalar field structure
  !> @param[in] u, v        - arrays containing x, y velocities
  !> @param[in] cell_mesh   - Mesh structure
  !> @param[in] n_int_cells - number of cells in the interior of the mesh
  !> @param[in,out] M       - Matrix structure being assigned
  subroutine compute_interior_coeffs(phi, u, v, cell_mesh, n_int_cells, M)
    use constants, only : add_mode
    use types, only: matrix_values, cell_locator, face_locator, neighbour_locator
    use utils, only: pack_entries, set_values
    use meshing, only: set_cell_location, set_face_location, set_neighbour_location, &
                       get_global_index, count_neighbours, get_boundary_status

    class(field), intent(in) :: phi
    real(accs_real), dimension(:), intent(in) :: u, v
    type(mesh), intent(in) :: cell_mesh
    integer(accs_int), intent(in) :: n_int_cells
    class(matrix), allocatable :: M

    type(matrix_values) :: mat_coeffs
    type(cell_locator) :: self_loc
    type(neighbour_locator) :: ngb_loc
    type(face_locator) :: face_loc
    integer(accs_int) :: self_idx, ngb_idx, local_idx
    integer(accs_int) :: j
    integer(accs_int) :: mat_counter
    integer(accs_int) :: n_ngb
    real(accs_real) :: face_area
    real(accs_real) :: diff_coeff, diff_coeff_total
    real(accs_real) :: adv_coeff, adv_coeff_total
    real(accs_real), dimension(ndim) :: face_normal
    logical :: is_boundary

    real(accs_real) :: mf !> The mass flux at a face

    mat_coeffs%mode = add_mode

    allocate(mat_coeffs%rglob(1))
    allocate(mat_coeffs%cglob(n_int_cells))
    allocate(mat_coeffs%val(n_int_cells))

    do local_idx = 1, cell_mesh%nlocal
      ! Calculate contribution from neighbours
      call set_cell_location(self_loc, cell_mesh, local_idx)
      call get_global_index(self_loc, self_idx)
      call count_neighbours(self_loc, n_ngb)
      mat_counter = 1
      adv_coeff_total = 0.0_accs_real
      diff_coeff_total = 0.0_accs_real
      do j = 1, n_ngb
        call set_neighbour_location(ngb_loc, self_loc, j)
        call get_global_index(ngb_loc, ngb_idx)
        call get_boundary_status(ngb_loc, is_boundary)

        diff_coeff = calc_diffusion_coeff(local_idx, j, cell_mesh)
        if (.not. is_boundary) then
          call set_face_location(face_loc, cell_mesh, local_idx, j)
          call get_face_area(face_loc, face_area)
          call get_face_normal(face_loc, face_normal)

          ! XXX: Why won't Fortran interfaces distinguish on extended types...
          ! TODO: This will be expensive (in a tight loop) - investigate moving to a type-bound
          !       procedure (should also eliminate the type check).
          mf = calc_mass_flux(u, v, ngb_idx, self_idx, face_area, face_normal, 0)
          select type(phi)
            type is(central_field)
              call calc_advection_coeff(phi, mf, 0, adv_coeff)
            type is(upwind_field)
              call calc_advection_coeff(phi, mf, 0, adv_coeff)
            class default
              print *, 'invalid velocity field discretisation'
              stop
          end select
          call pack_entries(mat_coeffs, 1, mat_counter, self_idx, ngb_idx, adv_coeff * mf + diff_coeff)
          mat_counter = mat_counter + 1
          adv_coeff_total = adv_coeff_total + (1.0_accs_real - adv_coeff) * mf
          diff_coeff_total = diff_coeff_total + diff_coeff
        else
          call pack_entries(mat_coeffs, 1, mat_counter, self_idx, -1, 0.0_accs_real)
          mat_counter = mat_counter + 1
        end if
      end do
      call pack_entries(mat_coeffs, 1, mat_counter, self_idx, self_idx, adv_coeff_total - diff_coeff_total)
      mat_counter = mat_counter + 1
      call set_values(mat_coeffs, M)
    end do

    deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)
  end subroutine compute_interior_coeffs

  !> @brief Computes the value of the scalar field on the boundary based on linear interpolation between 
  !  values provided on box corners
  !
  !> @param[in] ngb_index - index of neighbour with respect to CV (i.e. range 1-4 in square mesh)
  !> @param[in] row       - global row of cell within square mesh
  !> @param[in] col       - global column of cell within square mesh
  !> @param[in] cps       - number of cells per side in square mesh
  !> @param[in] bcs       - BC configuration data structure
  !> @param[out] bc_value - the value of the scalar field at the specified boundary
  subroutine compute_boundary_values(ngb_index, row, col, cps, bcs, bc_value)
    integer, intent(in) :: ngb_index  ! This is the index wrt the CV, not the ngb's cell index (i.e. range 1-4 for a square mesh)
    integer, intent(in) :: row
    integer, intent(in) :: col
    integer, intent(in) :: cps
    type(bc_config), intent(in) :: bcs
    real(accs_real), intent(out) :: bc_value
    real(accs_real) :: row_cps, col_cps

    row_cps = real(row, accs_real)/real(cps, accs_real)
    col_cps = real(col, accs_real)/real(cps, accs_real)

    bc_value = 0.0_accs_real
    if (bcs%bc_type(ngb_index) == bc_type_dirichlet .and. &
       (bcs%region(ngb_index) == bc_region_left .or. &
       bcs%region(ngb_index) == bc_region_right)) then
      bc_value = -((1.0_accs_real - row_cps) * bcs%endpoints(ngb_index, 1) + row_cps * bcs%endpoints(ngb_index, 2))
    else if (bcs%bc_type(ngb_index) == bc_type_dirichlet .and. &
            (bcs%region(ngb_index) == bc_region_top .or. &
            bcs%region(ngb_index) == bc_region_bottom)) then
      bc_value = -((1.0_accs_real - col_cps) * bcs%endpoints(ngb_index, 1) + col_cps * bcs%endpoints(ngb_index, 2))
    end if
  end subroutine compute_boundary_values

  !> @brief Computes the matrix coefficient for cells on the boundary of the mesh
  !
  !> @param[in] phi         - scalar field structure
  !> @param[in] u, v        - arrays containing x, y velocities
  !> @param[in] cell_mesh   - Mesh structure
  !> @param[in] bcs         - boundary conditions structure
  !> @param[in] cps         - number of cells per side
  !> @param[in,out] M       - Matrix structure being assigned
  !> @param[in,out] b       - vector structure being assigned
  subroutine compute_boundary_coeffs(phi, u, v, cell_mesh, bcs, cps, M, b)
    use constants, only : insert_mode, add_mode
    use types, only: matrix_values, vector_values, cell_locator, face_locator, neighbour_locator
    use utils, only: pack_entries, set_values
    use meshing, only: get_global_index, count_neighbours, get_boundary_status, &
                       set_cell_location, set_face_location, set_neighbour_location
    use bc_constants

    class(field), intent(in) :: phi
    real(accs_real), dimension(:), intent(in) :: u, v
    type(mesh), intent(in) :: cell_mesh
    type(bc_config), intent(in) :: bcs
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
    real(accs_real) :: face_area
    real(accs_real) :: diff_coeff
    real(accs_real) :: adv_coeff
    real(accs_real) :: bc_value
    real(accs_real), dimension(ndim) :: face_normal
    logical :: is_boundary

    real(accs_real) :: mf !> The mass flux at a face
    
    mat_coeffs%mode = add_mode
    b_coeffs%mode = add_mode

    allocate(mat_coeffs%rglob(1))
    allocate(mat_coeffs%cglob(1))
    allocate(mat_coeffs%val(1))
    allocate(b_coeffs%idx(1))
    allocate(b_coeffs%val(1))

    bc_counter = 1
    do local_idx = 1, cell_mesh%nlocal
      call set_cell_location(self_loc, cell_mesh, local_idx)
      call get_global_index(self_loc, self_idx)
      call count_neighbours(self_loc, n_ngb)
      ! Calculate contribution from neighbours
      do j = 1, n_ngb
        call set_neighbour_location(ngb_loc, self_loc, j)
        call get_global_index(ngb_loc, ngb_idx)
        call get_boundary_status(ngb_loc, is_boundary)

        call set_face_location(face_loc, cell_mesh, local_idx, j)
        call get_face_area(face_loc, face_area)
        call get_face_normal(face_loc, face_normal)

        mesh_ngb_idx = cell_mesh%nbidx(j, local_idx)
        if (is_boundary) then
          diff_coeff = calc_diffusion_coeff(local_idx, j, cell_mesh)
          mf = calc_mass_flux(u, v, ngb_idx, self_idx, face_area, face_normal, mesh_ngb_idx)
          select type(phi)
            type is(central_field)
              call calc_advection_coeff(phi, mf, mesh_ngb_idx, adv_coeff)
            type is(upwind_field)
              call calc_advection_coeff(phi, mf, mesh_ngb_idx, adv_coeff)
            class default
              print *, 'invalid velocity field discretisation'
              stop
          end select
          adv_coeff = adv_coeff * mf
          
          call calc_cell_coords(self_idx, cps, row, col)
          call compute_boundary_values(j, row, col, cps, bcs, bc_value)
          call pack_entries(b_coeffs, 1, self_idx, (adv_coeff + diff_coeff)*bc_value)
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
  !> @param[in] local_self_idx - the local cell index
  !> @param[in] local_ngb_idx  - the local neigbouring cell index
  !> @param[in] cell_mesh      - the mesh structure
  !> @param[out] coeff         - the diffusion coefficient
  module function calc_diffusion_coeff(local_self_idx, local_ngb_idx, cell_mesh) result(coeff)
    integer(accs_int), intent(in) :: local_self_idx
    integer(accs_int), intent(in) :: local_ngb_idx
    type(mesh), intent(in) :: cell_mesh
    real(accs_real) :: coeff

    type(face_locator) :: face_location
    real(accs_real) :: face_area
    real(accs_real), parameter :: diffusion_factor = 1.e-2_accs_real ! ALEXEI: temporarily hard-coded

    call set_face_location(face_location, cell_mesh, local_self_idx, local_ngb_idx)
    call get_face_area(face_location, face_area)

    coeff = -face_area*diffusion_factor/cell_mesh%h
  end function calc_diffusion_coeff

  !> @brief Calculates mass flux across given face. Note: assumes rho = 1 and uniform grid
  !
  !> @param[in] u, v     - arrays containing x, y velocities
  !> @param[in] ngb_idx  - Row and column of given neighbour in mesh
  !> @param[in] self_idx - Row and column of given cell in mesh
  !> @param[in] bc_flag  - Flag to indicate if neighbour is a boundary cell
  !> @param[out] flux    - The flux across the boundary
  module function calc_mass_flux(u, v, ngb_idx, self_idx, face_area, face_normal, bc_flag) result(flux)
    real(accs_real), dimension(:), intent(in) :: u, v
    integer(accs_int), intent(in) :: ngb_idx
    integer(accs_int), intent(in) :: self_idx
    real(accs_real), intent(in) :: face_area
    real(accs_real), dimension(ndim), intent(in) :: face_normal
    integer(accs_int), intent(in) :: bc_flag

    real(accs_real) :: flux

    flux = 0.0_accs_real
    
    ! ALEXEI: Write more general implementation handling BCs
    if (bc_flag == 0) then
      flux = 0.5_accs_real*(u(ngb_idx) + u(self_idx)) * face_area * face_normal(1) + &
             0.5_accs_real*(v(ngb_idx) + v(self_idx)) * face_area * face_normal(2)
    else if (bc_flag == -1 .or. bc_flag == -2) then
      flux = u(self_idx) * face_area
    else if (bc_flag == -3 .or. bc_flag == -4) then
      flux = v(self_idx) * face_area
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
