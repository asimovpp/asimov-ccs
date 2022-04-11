!> @brief Submodule file fv_CDS.smod
!> @build CDS
!
!> @details An implementation of the finite volume method using the CDS scheme

submodule (fv) fv_common

  use constants, only: add_mode, insert_mode, ndim
  use types, only: vector_values, matrix_values, cell_locator, face_locator, &
                   neighbour_locator
  use vec, only: get_vector_data, restore_vector_data
  use utils, only: pack_entries, set_values, update
  use meshing, only: count_neighbours, get_boundary_status, set_neighbour_location, &
                      get_local_index, get_global_index, get_volume, get_distance, &
                      set_face_location, get_face_area, get_face_normal, set_cell_location

  implicit none

contains

  !> @brief Computes fluxes and assign to matrix and RHS
  !
  !> @param[in] phi       - scalar field structure
  !> @param[in] mf        - mass flux field structure
  !> @param[in] cell_mesh - the mesh being used
  !> @param[in] bcs       - the boundary conditions structure being used
  !> @param[in] cps       - the number of cells per side in the (square) mesh
  !> @param[in,out] M     - Data structure containing matrix to be filled
  !> @param[in,out] vec   - Data structure containing RHS vector to be filled
  module subroutine compute_fluxes(phi, mf, cell_mesh, bcs, cps, M, vec)
    class(field), intent(in) :: phi
    class(field), intent(in) :: mf
    type(mesh), intent(in) :: cell_mesh
    type(bc_config), intent(in) :: bcs
    integer(accs_int), intent(in) :: cps
    class(matrix), intent(inout) :: M   
    class(vector), intent(inout) :: vec   

    integer(accs_int) :: n_int_cells
    real(accs_real), dimension(:), pointer :: mf_data

    associate (mf_vec => mf%vec)
      print *, "CF: get mf"
      call get_vector_data(mf_vec, mf_data)

      ! Loop over cells computing advection and diffusion fluxes
      n_int_cells = calc_matrix_nnz()
      print *, "CF: interior"
      call compute_interior_coeffs(phi, mf_data, cell_mesh, n_int_cells, M)

      ! Loop over boundaries
      print *, "CF: boundaries"
      call compute_boundary_coeffs(phi, mf_data, cell_mesh, bcs, cps, M, vec)

      print *, "CF: restore mf"
      call restore_vector_data(mf_vec, mf_data)
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
  !> @param[in] mf          - mass flux array defined at faces
  !> @param[in] cell_mesh   - Mesh structure
  !> @param[in] n_int_cells - number of cells in the interior of the mesh
  !> @param[in,out] M       - Matrix structure being assigned
  subroutine compute_interior_coeffs(phi, mf, cell_mesh, n_int_cells, M)
    class(field), intent(in) :: phi
    real(accs_real), dimension(:), intent(in) :: mf
    type(mesh), intent(in) :: cell_mesh
    integer(accs_int), intent(in) :: n_int_cells
    class(matrix), intent(inout) :: M

    type(matrix_values) :: mat_coeffs
    type(cell_locator) :: self_loc
    type(neighbour_locator) :: ngb_loc
    type(face_locator) :: face_loc
    integer(accs_int) :: self_idx, ngb_idx, local_idx, ngb_local_idx
    integer(accs_int) :: j
    integer(accs_int) :: mat_counter
    integer(accs_int) :: n_ngb
    real(accs_real) :: face_area
    real(accs_real) :: diff_coeff, diff_coeff_total
    real(accs_real) :: adv_coeff, adv_coeff_total
    real(accs_real), dimension(ndim) :: face_normal
    logical :: is_boundary

    integer(accs_int) :: idxf

    real(accs_real) :: sgn !> Sign indicating face orientation

    mat_coeffs%mode = add_mode

    allocate(mat_coeffs%rglob(1))
    allocate(mat_coeffs%cglob(n_int_cells))
    allocate(mat_coeffs%val(n_int_cells))

    do local_idx = 1, cell_mesh%nlocal
      ! Calculate contribution from neighbours
      call set_cell_location(cell_mesh, local_idx, self_loc)
      call get_global_index(self_loc, self_idx)
      call count_neighbours(self_loc, n_ngb)
      mat_counter = 1
      adv_coeff_total = 0.0_accs_real
      diff_coeff_total = 0.0_accs_real
      do j = 1, n_ngb
        call set_neighbour_location(self_loc, j, ngb_loc)
        call get_boundary_status(ngb_loc, is_boundary)

        if (.not. is_boundary) then
          diff_coeff = calc_diffusion_coeff(local_idx, j, cell_mesh)

          call get_global_index(ngb_loc, ngb_idx)
          call get_local_index(ngb_loc, ngb_local_idx)

          call set_face_location(cell_mesh, local_idx, j, face_loc)
          call get_face_area(face_loc, face_area)
          call get_face_normal(face_loc, face_normal)
          call get_local_index(face_loc, idxf)

          ! XXX: Why won't Fortran interfaces distinguish on extended types...
          ! TODO: This will be expensive (in a tight loop) - investigate moving to a type-bound
          !       procedure (should also eliminate the type check).
          if (ngb_local_idx < local_idx) then
            sgn = -1.0_accs_real
          else
            sgn = 1.0_accs_real
          end if
          select type(phi)
            type is(central_field)
              call calc_advection_coeff(phi, sgn * mf(idxf), 0, adv_coeff)
            type is(upwind_field)
              call calc_advection_coeff(phi, sgn * mf(idxf), 0, adv_coeff)
            class default
              print *, 'invalid velocity field discretisation'
              stop
          end select

          ! XXX: we are relying on div(u)=0 => a_P = -sum_nb a_nb
          adv_coeff = adv_coeff * (sgn * mf(idxf) * face_area)
          
          call pack_entries(1, mat_counter, self_idx, ngb_idx, adv_coeff + diff_coeff, mat_coeffs)
          mat_counter = mat_counter + 1
          adv_coeff_total = adv_coeff_total + adv_coeff
          diff_coeff_total = diff_coeff_total + diff_coeff
        else
          call pack_entries(1, mat_counter, self_idx, -1, 0.0_accs_real, mat_coeffs)
          mat_counter = mat_counter + 1
        end if
      end do
      call pack_entries(1, mat_counter, self_idx, self_idx, -(adv_coeff_total + diff_coeff_total), mat_coeffs)
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
    ! if (bcs%bc_type(ngb_index) == bc_type_dirichlet .and. &
    !    (bcs%region(ngb_index) == bc_region_left .or. &
    !    bcs%region(ngb_index) == bc_region_right)) then
    !   bc_value = -((1.0_accs_real - row_cps) * bcs%endpoints(ngb_index, 1) + row_cps * bcs%endpoints(ngb_index, 2))
    ! else if (bcs%bc_type(ngb_index) == bc_type_dirichlet .and. &
    !         (bcs%region(ngb_index) == bc_region_top .or. &
    !         bcs%region(ngb_index) == bc_region_bottom)) then
    !   bc_value = -((1.0_accs_real - col_cps) * bcs%endpoints(ngb_index, 1) + col_cps * bcs%endpoints(ngb_index, 2))
    ! end if

    if (bcs%bc_type(ngb_index) == 0) then
      bc_value = 0.0_accs_real
    else if (bcs%bc_type(ngb_index) == 1) then
      bc_value = 1.0_accs_real ! XXX: might not be correct
    else
      print *, "ERROR: Unknown boundary type ", bcs%bc_type(ngb_index)
    end if
    
  end subroutine compute_boundary_values

  !> @brief Computes the matrix coefficient for cells on the boundary of the mesh
  !
  !> @param[in] phi         - scalar field structure
  !> @param[in] mf          - mass flux array defined at faces
  !> @param[in] cell_mesh   - Mesh structure
  !> @param[in] bcs         - boundary conditions structure
  !> @param[in] cps         - number of cells per side
  !> @param[in,out] M       - Matrix structure being assigned
  !> @param[in,out] b       - vector structure being assigned
  subroutine compute_boundary_coeffs(phi, mf, cell_mesh, bcs, cps, M, b)
    class(field), intent(in) :: phi
    real(accs_real), dimension(:), intent(in) :: mf
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
    integer(accs_int) :: self_idx, local_idx
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

    integer(accs_int) :: idxf
    
    mat_coeffs%mode = add_mode
    b_coeffs%mode = add_mode

    allocate(mat_coeffs%rglob(1))
    allocate(mat_coeffs%cglob(1))
    allocate(mat_coeffs%val(1))
    allocate(b_coeffs%idx(1))
    allocate(b_coeffs%val(1))

    bc_counter = 1
    do local_idx = 1, cell_mesh%nlocal
      call set_cell_location(cell_mesh, local_idx, self_loc)
      call get_global_index(self_loc, self_idx)
      call count_neighbours(self_loc, n_ngb)
      ! Calculate contribution from neighbours
      do j = 1, n_ngb
        call set_neighbour_location(self_loc, j, ngb_loc)
        call get_boundary_status(ngb_loc, is_boundary)
        if (is_boundary) then
          ! call get_global_index(ngb_loc, ngb_idx)
          call get_local_index(ngb_loc, mesh_ngb_idx)

          call set_face_location(cell_mesh, local_idx, j, face_loc)
          call get_face_area(face_loc, face_area)
          call get_face_normal(face_loc, face_normal)
          call get_local_index(face_loc, idxf)
          
          diff_coeff = calc_diffusion_coeff(local_idx, j, cell_mesh)
          select type(phi)
            type is(central_field)
              call calc_advection_coeff(phi, mf(idxf), mesh_ngb_idx, adv_coeff)
            type is(upwind_field)
              call calc_advection_coeff(phi, mf(idxf), mesh_ngb_idx, adv_coeff)
            class default
              print *, 'invalid velocity field discretisation'
              stop
          end select
          adv_coeff = adv_coeff * (mf(idxf) * face_area)

          call calc_cell_coords(self_idx, cps, row, col)
          call compute_boundary_values(j, row, col, cps, bcs, bc_value)
          call pack_entries(1, self_idx, -(adv_coeff + diff_coeff)*bc_value, b_coeffs)
          call pack_entries(1, 1, self_idx, self_idx, -(adv_coeff + diff_coeff), mat_coeffs)
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
    real(accs_real), parameter :: diffusion_factor = 1.e-2_accs_real ! XXX: temporarily hard-coded
    logical :: is_boundary
    real(accs_real), dimension(ndim) :: dx
    real(accs_real) :: dxmag
    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb

    call set_face_location(cell_mesh, local_self_idx, local_ngb_idx, face_location)
    call get_face_area(face_location, face_area)
    call get_boundary_status(face_location, is_boundary)

    call set_cell_location(cell_mesh, local_self_idx, loc_p)
    if (.not. is_boundary) then
      call set_neighbour_location(loc_p, local_ngb_idx, loc_nb)
      call get_distance(loc_p, loc_nb, dx)
    else
      call get_distance(loc_p, face_location, dx)
    end if
    dxmag = sqrt(sum(dx**2))
    
    coeff = -face_area * diffusion_factor / dxmag
  end function calc_diffusion_coeff

  !> @brief Calculates mass flux across given face. Note: assumes rho = 1 and uniform grid
  !
  !> @param[in] u, v     - arrays containing x, y velocities
  !> @param[in] p        - array containing pressure
  !> @param[in] pgradx   - array containing pressure gradient in x
  !> @param[in] pgrady   - array containing pressure gradient in y
  !> @param[in] invAu    - array containing inverse momentum diagonal in x
  !> @param[in] invAv    - array containing inverse momentum diagonal in y
  !> @param[in] loc_f    - face locator
  !> @param[out] flux    - The flux across the boundary
  module function calc_mass_flux(u, v, p, pgradx, pgrady, invAu, invAv, loc_f) result(flux)
    real(accs_real), dimension(:), intent(in) :: u, v
    real(accs_real), dimension(:), intent(in) :: p
    real(accs_real), dimension(:), intent(in) :: pgradx, pgrady
    real(accs_real), dimension(:), intent(in) :: invAu, invAv
    type(face_locator), intent(in) :: loc_f

    real(accs_real) :: flux

    ! Local variables
    logical :: is_boundary                          !> Boundary indicator
    type(cell_locator) :: loc_p                     !> Primary cell locator
    type(neighbour_locator) :: loc_nb               !> Neighbour cell locator
    integer(accs_int) :: idxnb                      !> Neighbour cell index
    real(accs_real) :: flux_corr                    !> Flux correction
    real(accs_real), dimension(ndim) :: dx          !> Cell-cell distance
    real(accs_real) :: dxmag                        !> Cell-cell distance magnitude
    real(accs_real), dimension(ndim) :: face_normal !> (local) face-normal array
    real(accs_real) :: Vp                           !> Primary cell volume
    real(accs_real) :: Vnb                          !> Neighbour cell volume
    real(accs_real) :: Vf                           !> Face "volume"
    real(accs_real) :: invAp                        !> Primary cell inverse momentum coefficient
    real(accs_real) :: invAnb                       !> Neighbour cell inverse momentum coefficient
    real(accs_real) :: invAf                        !> Face inverse momentum coefficient
    
    call get_boundary_status(loc_f, is_boundary)
    if (.not. is_boundary) then
      associate(mesh => loc_f%mesh, &
           idxp => loc_f%cell_idx, &
           j => loc_f%cell_face_ctr)
        
        call set_cell_location(mesh, idxp, loc_p)
        call set_neighbour_location(loc_p, j, loc_nb)
        call get_local_index(loc_nb, idxnb)

        call get_face_normal(loc_f, face_normal)
        
        flux = 0.5_accs_real * ((u(idxp) + u(idxnb)) * face_normal(1) &
             + (v(idxp) + v(idxnb)) * face_normal(2))

        !
        ! Rhie-Chow correction from Ferziger & Peric
        !
        call get_distance(loc_p, loc_nb, dx)
        dxmag = sqrt(sum(dx**2))
        call get_face_normal(loc_f, face_normal)
        flux_corr = -(p(idxnb) - p(idxp)) / dxmag
        flux_corr = flux_corr + 0.5_accs_real * ((pgradx(idxp) + pgradx(idxnb)) * face_normal(1) &
             + (pgrady(idxp) + pgrady(idxnb)) * face_normal(2))

        call get_volume(loc_p, Vp)
        call get_volume(loc_nb, Vnb)
        Vf = 0.5_accs_real * (Vp + Vnb)

        ! This is probably not quite right ...
        invAp = 0.5_accs_real * (invAu(idxp) + invAv(idxp))
        invAnb = 0.5_accs_real * (invAu(idxnb) + invAv(idxnb))
        invAf = 0.5_accs_real * (invAp + invAnb)
        
        flux_corr = (Vf * invAf) * flux_corr
          
        ! Apply correction
        flux = flux + flux_corr

        if (idxp > idxnb) then
          ! XXX: making convention to point from low to high cell!
          flux = -flux
        end if
      end associate
    else 
      ! TODO: Write more general implementation handling BCs
      flux = 0.0_accs_real ! XXX: hardcoded zero-flux BC
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

  !> @brief Performs an update of the gradients of a field.
  !
  !> @param[in]    cell_mesh - the mesh
  !> @param[inout] phi       - the field whose gradients we want to update
  !
  !> @note This will perform a parallel update of the gradient fields to ensure halo cells are
  !!       correctly updated on other PEs.
  module subroutine update_gradient(cell_mesh, phi)

    type(mesh), intent(in) :: cell_mesh
    class(field), intent(inout) :: phi

    real(accs_real), dimension(:), pointer :: gradx_data, grady_data, gradz_data
    real(accs_real), dimension(:), allocatable :: gradx_old, grady_old, gradz_old

    integer(accs_real) :: i
    
    call get_vector_data(phi%gradx, gradx_data)
    call get_vector_data(phi%grady, grady_data)
    call get_vector_data(phi%gradz, gradz_data)

    associate(ntotal => cell_mesh%ntotal)
      allocate(gradx_old(ntotal))
      allocate(grady_old(ntotal))
      allocate(gradz_old(ntotal))
      do i = 1, ntotal
        gradx_old(i) = gradx_data(i)
        grady_old(i) = grady_data(i)
        gradz_old(i) = gradz_data(i)
      end do
    end associate
    
    call restore_vector_data(phi%gradx, gradx_data)
    call restore_vector_data(phi%grady, grady_data)
    call restore_vector_data(phi%gradz, gradz_data)
    
    call update_gradient_component(cell_mesh, 1, phi%vec, gradx_old, grady_old, gradz_old, phi%gradx)
    call update(phi%gradx) ! XXX: opportunity to overlap update with later compute (begin/compute/end)
    call update_gradient_component(cell_mesh, 2, phi%vec, gradx_old, grady_old, gradz_old, phi%grady)
    call update(phi%grady) ! XXX: opportunity to overlap update with later compute (begin/compute/end)
    call update_gradient_component(cell_mesh, 3, phi%vec, gradx_old, grady_old, gradz_old, phi%gradz)
    call update(phi%gradz) ! XXX: opportunity to overlap update with later compute (begin/compute/end)

    deallocate(gradx_old)
    deallocate(grady_old)
    deallocate(gradz_old)
    
  end subroutine update_gradient

  !> @brief Helper subroutine to calculate a gradient component at a time.
  !
  !> @param[in] cell_mesh   - the mesh
  !> @param[in] component   - which vector component (i.e. direction) to update?
  !> @param[in] phi         - a cell-centred array of the field whose gradient we
  !!                          want to compute
  !> @param[inout] gradient - a cell-centred array of the gradient
  subroutine update_gradient_component(cell_mesh, component, phi, gradx_old, grady_old, gradz_old, gradient)


    type(mesh), intent(in) :: cell_mesh
    integer(accs_int), intent(in) :: component
    class(vector), intent(in) :: phi
    real(accs_real), dimension(:), intent(in) :: gradx_old
    real(accs_real), dimension(:), intent(in) :: grady_old
    real(accs_real), dimension(:), intent(in) :: gradz_old
    class(vector), intent(inout) :: gradient
    
    type(vector_values) :: grad_values
    real(accs_real), dimension(:), pointer :: phi_data
    real(accs_real) :: grad
    
    integer(accs_int) :: i
    integer(accs_int) :: j
    type(cell_locator) :: loc_p
    type(face_locator) :: loc_f
    type(neighbour_locator) :: loc_nb
    
    integer(accs_int) :: nnb
    integer(accs_int) :: nb
    
    real(accs_real) :: phif

    logical :: is_boundary

    real(accs_real) :: face_area
    real(accs_real), dimension(ndim) :: face_norm

    real(accs_real) :: V
    integer(accs_int) :: idxg

    real(accs_real), dimension(ndim) :: dx
    
    allocate(grad_values%idx(1))
    allocate(grad_values%val(1))
    grad_values%mode = insert_mode

    call get_vector_data(phi, phi_data)
    
    do i = 1, cell_mesh%nlocal
      grad = 0.0_accs_int
      
      call set_cell_location(cell_mesh, i, loc_p)
      call count_neighbours(loc_p, nnb)
      do j = 1, nnb
        call set_face_location(cell_mesh, i, j, loc_f)
        call get_boundary_status(loc_f, is_boundary)

        if (.not. is_boundary) then
          call set_neighbour_location(loc_p, j, loc_nb)
          call get_local_index(loc_nb, nb)
          phif = 0.5_accs_real * (phi_data(i) + phi_data(nb)) ! XXX: Need to do proper interpolation
        else
          call get_distance(loc_p, loc_f, dx)
          phif = phi_data(i) + (gradx_old(i) * dx(1) + grady_old(i) * dx(2) + gradz_old(i) * dx(3))
        end if

        call get_face_area(loc_f, face_area)
        call get_face_normal(loc_f, face_norm)

        grad = grad + phif * (face_area * face_norm(component))
      end do

      call get_volume(loc_p, V)
      grad = grad / V
      
      call get_global_index(loc_p, idxg)
      call pack_entries(1, idxg, grad, grad_values)
      call set_values(grad_values, gradient)
    end do

    call restore_vector_data(phi, phi_data)
    
  end subroutine update_gradient_component
  
end submodule fv_common
