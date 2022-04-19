!>  Submodule file fv_CDS.smod
!>  CDS
!
!>  An implementation of the finite volume method using the CDS scheme

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

  !>  Computes fluxes and assign to matrix and RHS
  module subroutine compute_fluxes(phi, mf, mesh, bcs, cps, M, vec)
    class(field), intent(in) :: phi           !< scalar field structure
    class(field), intent(in) :: mf            !< mass flux field structure
    type(ccs_mesh), intent(in) :: mesh        !< the mesh being used
    type(bc_config), intent(in) :: bcs        !< the boundary conditions structure being used
    integer(ccs_int), intent(in) :: cps       !< the number of cells per side in the (square) mesh
    class(ccs_matrix), intent(inout) :: M     !< Data structure containing matrix to be filled
    class(ccs_vector), intent(inout) :: vec   !< Data structure containing RHS vector to be filled

    integer(ccs_int) :: n_int_cells
    real(ccs_real), dimension(:), pointer :: mf_data

    associate (mf_values => mf%values)
      print *, "CF: get mf"
      call get_vector_data(mf_values, mf_data)

      ! Loop over cells computing advection and diffusion fluxes
      n_int_cells = calc_matrix_nnz()
      print *, "CF: interior"
      call compute_interior_coeffs(phi, mf_data, mesh, n_int_cells, M)

      ! Loop over boundaries
      print *, "CF: boundaries"
      call compute_boundary_coeffs(phi, mf_data, mesh, bcs, cps, M, vec)

      print *, "CF: restore mf"
      call restore_vector_data(mf_values, mf_data)
    end associate

  end subroutine compute_fluxes

  !>  Returns the number of entries per row that are non-zero
  !
  !>  This assumes a square 2d grid
  pure function calc_matrix_nnz() result(nnz)
    integer(ccs_int) :: nnz   !< number of non-zero entries per row 

    nnz = 5_ccs_int
  end function calc_matrix_nnz

  !>  Computes the matrix coefficient for cells in the interior of the mesh
  subroutine compute_interior_coeffs(phi, mf, mesh, n_int_cells, M)
    class(field), intent(in) :: phi                 !< scalar field structure
    real(ccs_real), dimension(:), intent(in) :: mf  !< mass flux array defined at faces
    type(ccs_mesh), intent(in) :: mesh              !< Mesh structure
    integer(ccs_int), intent(in) :: n_int_cells     !< number of cells in the interior of the mesh
    class(ccs_matrix), intent(inout) :: M           !< Matrix structure being assigned

    type(matrix_values) :: mat_coeffs
    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb
    type(face_locator) :: loc_f
    integer(ccs_int) :: global_index_p, global_index_nb, index_p, index_nb
    integer(ccs_int) :: j
    integer(ccs_int) :: mat_counter
    integer(ccs_int) :: nnb
    real(ccs_real) :: face_area
    real(ccs_real) :: diff_coeff, diff_coeff_total
    real(ccs_real) :: adv_coeff, adv_coeff_total
    real(ccs_real), dimension(ndim) :: face_normal
    logical :: is_boundary

    integer(ccs_int) :: index_f

    real(ccs_real) :: sgn !< Sign indicating face orientation

    mat_coeffs%setter_mode = add_mode

    allocate(mat_coeffs%row_indices(1))
    allocate(mat_coeffs%col_indices(n_int_cells))
    allocate(mat_coeffs%values(n_int_cells))

    do index_p = 1, mesh%nlocal
      ! Calculate contribution from neighbours
      call set_cell_location(mesh, index_p, loc_p)
      call get_global_index(loc_p, global_index_p)
      call count_neighbours(loc_p, nnb)
      mat_counter = 1
      adv_coeff_total = 0.0_ccs_real
      diff_coeff_total = 0.0_ccs_real
      do j = 1, nnb
        call set_neighbour_location(loc_p, j, loc_nb)
        call get_boundary_status(loc_nb, is_boundary)

        if (.not. is_boundary) then
          diff_coeff = calc_diffusion_coeff(index_p, j, mesh)

          call get_global_index(loc_nb, global_index_nb)
          call get_local_index(loc_nb, index_nb)

          call set_face_location(mesh, index_p, j, loc_f)
          call get_face_area(loc_f, face_area)
          call get_face_normal(loc_f, face_normal)
          call get_local_index(loc_f, index_f)

          ! XXX: Why won't Fortran interfaces distinguish on extended types...
          ! TODO: This will be expensive (in a tight loop) - investigate moving to a type-bound
          !       procedure (should also eliminate the type check).
          if (index_nb < index_p) then
            sgn = -1.0_ccs_real
          else
            sgn = 1.0_ccs_real
          end if
          select type(phi)
            type is(central_field)
              call calc_advection_coeff(phi, sgn * mf(index_f), 0, adv_coeff)
            type is(upwind_field)
              call calc_advection_coeff(phi, sgn * mf(index_f), 0, adv_coeff)
            class default
              print *, 'invalid velocity field discretisation'
              stop
          end select

          ! XXX: we are relying on div(u)=0 => a_P = -sum_nb a_nb
          adv_coeff = adv_coeff * (sgn * mf(index_f) * face_area)
          
          call pack_entries(1, mat_counter, global_index_p, global_index_nb, adv_coeff + diff_coeff, mat_coeffs)
          mat_counter = mat_counter + 1
          adv_coeff_total = adv_coeff_total + adv_coeff
          diff_coeff_total = diff_coeff_total + diff_coeff
        else
          call pack_entries(1, mat_counter, global_index_p, -1, 0.0_ccs_real, mat_coeffs)
          mat_counter = mat_counter + 1
        end if
      end do
      call pack_entries(1, mat_counter, global_index_p, global_index_p, -(adv_coeff_total + diff_coeff_total), mat_coeffs)
      mat_counter = mat_counter + 1
      call set_values(mat_coeffs, M)
    end do

    deallocate(mat_coeffs%row_indices)
    deallocate(mat_coeffs%col_indices)
    deallocate(mat_coeffs%values)
  end subroutine compute_interior_coeffs

  !v  Computes the value of the scalar field on the boundary based on linear interpolation between 
  !  values provided on box corners
  subroutine compute_boundary_values(index_nb, row, col, cps, bcs, bc_value)
    integer, intent(in) :: index_nb           !< index of neighbour with respect to CV (i.e. range 1-4 in square mesh)
    integer, intent(in) :: row                !< global row of cell within square mesh
    integer, intent(in) :: col                !< global column of cell within square mesh
    integer, intent(in) :: cps                !< number of cells per side in square mesh
    type(bc_config), intent(in) :: bcs        !< BC configuration data structure
    real(ccs_real), intent(out) :: bc_value   !< the value of the scalar field at the specified boundary
    real(ccs_real) :: row_cps, col_cps

    row_cps = real(row, ccs_real)/real(cps, ccs_real)
    col_cps = real(col, ccs_real)/real(cps, ccs_real)

    bc_value = 0.0_ccs_real
    ! if (bcs%bc_type(index_nb) == bc_type_dirichlet .and. &
    !    (bcs%region(index_nb) == bc_region_left .or. &
    !    bcs%region(index_nb) == bc_region_right)) then
    !   bc_value = -((1.0_ccs_real - row_cps) * bcs%endpoints(index_nb, 1) + row_cps * bcs%endpoints(index_nb, 2))
    ! else if (bcs%bc_type(index_nb) == bc_type_dirichlet .and. &
    !         (bcs%region(index_nb) == bc_region_top .or. &
    !         bcs%region(index_nb) == bc_region_bottom)) then
    !   bc_value = -((1.0_ccs_real - col_cps) * bcs%endpoints(index_nb, 1) + col_cps * bcs%endpoints(index_nb, 2))
    ! end if

    if (bcs%bc_type(index_nb) == 0) then
      bc_value = 0.0_ccs_real
    else if (bcs%bc_type(index_nb) == 1) then
      bc_value = 1.0_ccs_real ! XXX: might not be correct
    else
      print *, "ERROR: Unknown boundary type ", bcs%bc_type(index_nb)
    end if
    
  end subroutine compute_boundary_values

  !>  Computes the matrix coefficient for cells on the boundary of the mesh
  subroutine compute_boundary_coeffs(phi, mf, mesh, bcs, cps, M, b)
    class(field), intent(in) :: phi                 !< scalar field structure
    real(ccs_real), dimension(:), intent(in) :: mf  !< mass flux array defined at faces
    type(ccs_mesh), intent(in) :: mesh              !< Mesh structure
    type(bc_config), intent(in) :: bcs              !< boundary conditions structure
    integer(ccs_int), intent(in) :: cps             !< number of cells per side
    class(ccs_matrix), intent(inout) :: M           !< Matrix structure being assigned
    class(ccs_vector), intent(inout) :: b           !< vector structure being assigned

    type(matrix_values) :: mat_coeffs
    type(vector_values) :: b_coeffs
    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb
    type(face_locator) :: loc_f
    integer(ccs_int) :: global_index_p, index_p
    integer(ccs_int) :: j
    integer(ccs_int) :: bc_counter
    integer(ccs_int) :: row, col
    integer(ccs_int) :: nnb, index_nb
    real(ccs_real) :: face_area
    real(ccs_real) :: diff_coeff
    real(ccs_real) :: adv_coeff
    real(ccs_real) :: bc_value
    real(ccs_real), dimension(ndim) :: face_normal
    logical :: is_boundary

    integer(ccs_int) :: index_f
    
    mat_coeffs%setter_mode = add_mode
    b_coeffs%setter_mode = add_mode

    allocate(mat_coeffs%row_indices(1))
    allocate(mat_coeffs%col_indices(1))
    allocate(mat_coeffs%values(1))
    allocate(b_coeffs%indices(1))
    allocate(b_coeffs%values(1))

    bc_counter = 1
    do index_p = 1, mesh%nlocal
      call set_cell_location(mesh, index_p, loc_p)
      call get_global_index(loc_p, global_index_p)
      call count_neighbours(loc_p, nnb)
      ! Calculate contribution from neighbours
      do j = 1, nnb
        call set_neighbour_location(loc_p, j, loc_nb)
        call get_boundary_status(loc_nb, is_boundary)
        if (is_boundary) then
          ! call get_global_index(loc_nb, global_index_nb)
          call get_local_index(loc_nb, index_nb)

          call set_face_location(mesh, index_p, j, loc_f)
          call get_face_area(loc_f, face_area)
          call get_face_normal(loc_f, face_normal)
          call get_local_index(loc_f, index_f)
          
          diff_coeff = calc_diffusion_coeff(index_p, j, mesh)
          select type(phi)
            type is(central_field)
              call calc_advection_coeff(phi, mf(index_f), index_nb, adv_coeff)
            type is(upwind_field)
              call calc_advection_coeff(phi, mf(index_f), index_nb, adv_coeff)
            class default
              print *, 'invalid velocity field discretisation'
              stop
          end select
          adv_coeff = adv_coeff * (mf(index_f) * face_area)

          call calc_cell_coords(global_index_p, cps, row, col)
          call compute_boundary_values(j, row, col, cps, bcs, bc_value)
          call pack_entries(1, global_index_p, -(adv_coeff + diff_coeff)*bc_value, b_coeffs)
          call pack_entries(1, 1, global_index_p, global_index_p, -(adv_coeff + diff_coeff), mat_coeffs)
          call set_values(b_coeffs, b)
          call set_values(mat_coeffs, M)
          bc_counter = bc_counter + 1
        end if
      end do
    end do
    deallocate(mat_coeffs%row_indices)
    deallocate(mat_coeffs%col_indices)
    deallocate(mat_coeffs%values)
    deallocate(b_coeffs%indices)
    deallocate(b_coeffs%values)

  end subroutine compute_boundary_coeffs

  !>  Sets the diffusion coefficient
  module function calc_diffusion_coeff(index_p, index_nb, mesh) result(coeff)
    integer(ccs_int), intent(in) :: index_p     !< the local cell index
    integer(ccs_int), intent(in) :: index_nb    !< the local neigbouring cell index
    type(ccs_mesh), intent(in) :: mesh          !< the mesh structure
    real(ccs_real) :: coeff                     !< the diffusion coefficient

    type(face_locator) :: loc_f
    real(ccs_real) :: face_area
    real(ccs_real), parameter :: diffusion_factor = 1.e-2_ccs_real ! XXX: temporarily hard-coded
    logical :: is_boundary
    real(ccs_real), dimension(ndim) :: dx
    real(ccs_real) :: dxmag
    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb

    call set_face_location(mesh, index_p, index_nb, loc_f)
    call get_face_area(loc_f, face_area)
    call get_boundary_status(loc_f, is_boundary)

    call set_cell_location(mesh, index_p, loc_p)
    if (.not. is_boundary) then
      call set_neighbour_location(loc_p, index_nb, loc_nb)
      call get_distance(loc_p, loc_nb, dx)
    else
      call get_distance(loc_p, loc_f, dx)
    end if
    dxmag = sqrt(sum(dx**2))
    
    coeff = -face_area * diffusion_factor / dxmag
  end function calc_diffusion_coeff

  !>  Calculates mass flux across given face. Note: assumes rho = 1 and uniform grid
  !
  !> @param[in] u, v     - arrays containing x, y velocities
  !> @param[in] p        - array containing pressure
  !> @param[in] dpdx   - array containing pressure gradient in x
  !> @param[in] dpdy   - array containing pressure gradient in y
  !> @param[in] invAu    - array containing inverse momentum diagonal in x
  !> @param[in] invAv    - array containing inverse momentum diagonal in y
  !> @param[in] loc_f    - face locator
  !> @param[out] flux    - The flux across the boundary
  module function calc_mass_flux(u, v, p, dpdx, dpdy, invAu, invAv, loc_f) result(flux)
    real(ccs_real), dimension(:), intent(in) :: u, v
    real(ccs_real), dimension(:), intent(in) :: p
    real(ccs_real), dimension(:), intent(in) :: dpdx, dpdy
    real(ccs_real), dimension(:), intent(in) :: invAu, invAv
    type(face_locator), intent(in) :: loc_f

    real(ccs_real) :: flux

    ! Local variables
    logical :: is_boundary                         !< Boundary indicator
    type(cell_locator) :: loc_p                    !< Primary cell locator
    type(neighbour_locator) :: loc_nb              !< Neighbour cell locator
    integer(ccs_int) :: index_nb                   !< Neighbour cell index
    real(ccs_real) :: flux_corr                    !< Flux correction
    real(ccs_real), dimension(ndim) :: dx          !< Cell-cell distance
    real(ccs_real) :: dxmag                        !< Cell-cell distance magnitude
    real(ccs_real), dimension(ndim) :: face_normal !< (local) face-normal array
    real(ccs_real) :: Vp                           !< Primary cell volume
    real(ccs_real) :: V_nb                         !< Neighbour cell volume
    real(ccs_real) :: Vf                           !< Face "volume"
    real(ccs_real) :: invAp                        !< Primary cell inverse momentum coefficient
    real(ccs_real) :: invA_nb                      !< Neighbour cell inverse momentum coefficient
    real(ccs_real) :: invAf                        !< Face inverse momentum coefficient
    
    call get_boundary_status(loc_f, is_boundary)
    if (.not. is_boundary) then
      associate(mesh => loc_f%mesh, &
           index_p => loc_f%index_p, &
           j => loc_f%cell_face_ctr)
        
        call set_cell_location(mesh, index_p, loc_p)
        call set_neighbour_location(loc_p, j, loc_nb)
        call get_local_index(loc_nb, index_nb)

        call get_face_normal(loc_f, face_normal)
        
        flux = 0.5_ccs_real * ((u(index_p) + u(index_nb)) * face_normal(1) &
             + (v(index_p) + v(index_nb)) * face_normal(2))

        !
        ! Rhie-Chow correction from Ferziger & Peric
        !
        call get_distance(loc_p, loc_nb, dx)
        dxmag = sqrt(sum(dx**2))
        call get_face_normal(loc_f, face_normal)
        flux_corr = -(p(index_nb) - p(index_p)) / dxmag
        flux_corr = flux_corr + 0.5_ccs_real * ((dpdx(index_p) + dpdx(index_nb)) * face_normal(1) &
             + (dpdy(index_p) + dpdy(index_nb)) * face_normal(2))

        call get_volume(loc_p, Vp)
        call get_volume(loc_nb, V_nb)
        Vf = 0.5_ccs_real * (Vp + V_nb)

        ! This is probably not quite right ...
        invAp = 0.5_ccs_real * (invAu(index_p) + invAv(index_p))
        invA_nb = 0.5_ccs_real * (invAu(index_nb) + invAv(index_nb))
        invAf = 0.5_ccs_real * (invAp + invA_nb)
        
        flux_corr = (Vf * invAf) * flux_corr
          
        ! Apply correction
        flux = flux + flux_corr

        if (index_p > index_nb) then
          ! XXX: making convention to point from low to high cell!
          flux = -flux
        end if
      end associate
    else 
      ! TODO: Write more general implementation handling BCs
      flux = 0.0_ccs_real ! XXX: hardcoded zero-flux BC
    end if
    
  end function calc_mass_flux

  !>  Calculates the row and column indices from flattened vector index. Assumes square mesh
  !
  !> @param[in] index  - cell index
  !> @param[in] cps  - number of cells per side
  !> @param[out] row - cell row within mesh
  !> @param[out] col - cell column within mesh
  module subroutine calc_cell_coords(index, cps, row, col)
    integer(ccs_int), intent(in) :: index, cps
    integer(ccs_int), intent(out) :: row, col

    col = modulo(index-1,cps) + 1 
    row = (index-1)/cps + 1
  end subroutine calc_cell_coords

  !>  Performs an update of the gradients of a field.
  !
  !> @param[in]    mesh - the mesh
  !> @param[inout] phi       - the field whose gradients we want to update
  !
  !> @note This will perform a parallel update of the gradient fields to ensure halo cells are
  !!       correctly updated on other PEs.

  module subroutine update_gradient(mesh, phi)
  
    type(ccs_mesh), intent(in) :: mesh
    class(field), intent(inout) :: phi

    real(ccs_real), dimension(:), pointer :: x_gradients_data, y_gradients_data, z_gradients_data
    real(ccs_real), dimension(:), allocatable :: x_gradients_old, y_gradients_old, z_gradients_old

    integer(ccs_real) :: i
    
    call get_vector_data(phi%x_gradients, x_gradients_data)
    call get_vector_data(phi%y_gradients, y_gradients_data)
    call get_vector_data(phi%z_gradients, z_gradients_data)

    associate(ntotal => mesh%ntotal)
      allocate(x_gradients_old(ntotal))
      allocate(y_gradients_old(ntotal))
      allocate(z_gradients_old(ntotal))
      do i = 1, ntotal
        x_gradients_old(i) = x_gradients_data(i)
        y_gradients_old(i) = y_gradients_data(i)
        z_gradients_old(i) = z_gradients_data(i)
      end do
    end associate
    
    call restore_vector_data(phi%x_gradients, x_gradients_data)
    call restore_vector_data(phi%y_gradients, y_gradients_data)
    call restore_vector_data(phi%z_gradients, z_gradients_data)
    
    call update_gradient_component(mesh, 1, phi%values, x_gradients_old, y_gradients_old, z_gradients_old, phi%x_gradients)
    call update(phi%x_gradients) ! XXX: opportunity to overlap update with later compute (begin/compute/end)
    call update_gradient_component(mesh, 2, phi%values, x_gradients_old, y_gradients_old, z_gradients_old, phi%y_gradients)
    call update(phi%y_gradients) ! XXX: opportunity to overlap update with later compute (begin/compute/end)
    call update_gradient_component(mesh, 3, phi%values, x_gradients_old, y_gradients_old, z_gradients_old, phi%z_gradients)
    call update(phi%z_gradients) ! XXX: opportunity to overlap update with later compute (begin/compute/end)

    deallocate(x_gradients_old)
    deallocate(y_gradients_old)
    deallocate(z_gradients_old)
    
  end subroutine update_gradient

  !>  Helper subroutine to calculate a gradient component at a time.
  !
  !> @param[in] mesh   - the mesh
  !> @param[in] component   - which vector component (i.e. direction) to update?
  !> @param[in] phi         - a cell-centred array of the field whose gradient we
  !!                          want to compute
  !> @param[inout] gradients - a cell-centred array of the gradient
  subroutine update_gradient_component(mesh, component, phi, x_gradients_old, y_gradients_old, z_gradients_old, gradients)


    type(ccs_mesh), intent(in) :: mesh
    integer(ccs_int), intent(in) :: component
    class(ccs_vector), intent(in) :: phi
    real(ccs_real), dimension(:), intent(in) :: x_gradients_old
    real(ccs_real), dimension(:), intent(in) :: y_gradients_old
    real(ccs_real), dimension(:), intent(in) :: z_gradients_old
    class(ccs_vector), intent(inout) :: gradients
    
    type(vector_values) :: grad_values
    real(ccs_real), dimension(:), pointer :: phi_data
    real(ccs_real) :: grad
    
    integer(ccs_int) :: i
    integer(ccs_int) :: j
    type(cell_locator) :: loc_p
    type(face_locator) :: loc_f
    type(neighbour_locator) :: loc_nb
    
    integer(ccs_int) :: nnb
    integer(ccs_int) :: nb
    
    real(ccs_real) :: phif

    logical :: is_boundary

    real(ccs_real) :: face_area
    real(ccs_real), dimension(ndim) :: face_norm

    real(ccs_real) :: V
    integer(ccs_int) :: global_index_p

    real(ccs_real), dimension(ndim) :: dx
    
    allocate(grad_values%indices(1))
    allocate(grad_values%values(1))
    grad_values%setter_mode = insert_mode

    call get_vector_data(phi, phi_data)
    
    do i = 1, mesh%nlocal
      grad = 0.0_ccs_int
      
      call set_cell_location(mesh, i, loc_p)
      call count_neighbours(loc_p, nnb)
      do j = 1, nnb
        call set_face_location(mesh, i, j, loc_f)
        call get_boundary_status(loc_f, is_boundary)

        if (.not. is_boundary) then
          call set_neighbour_location(loc_p, j, loc_nb)
          call get_local_index(loc_nb, nb)
          phif = 0.5_ccs_real * (phi_data(i) + phi_data(nb)) ! XXX: Need to do proper interpolation
        else
          call get_distance(loc_p, loc_f, dx)
          phif = phi_data(i) + (x_gradients_old(i) * dx(1) + y_gradients_old(i) * dx(2) + z_gradients_old(i) * dx(3))
        end if

        call get_face_area(loc_f, face_area)
        call get_face_normal(loc_f, face_norm)

        grad = grad + phif * (face_area * face_norm(component))
      end do

      call get_volume(loc_p, V)
      grad = grad / V
      
      call get_global_index(loc_p, global_index_p)
      call pack_entries(1, global_index_p, grad, grad_values)
      call set_values(grad_values, gradients)
    end do

    call restore_vector_data(phi, phi_data)
    
  end subroutine update_gradient_component
  
end submodule fv_common
