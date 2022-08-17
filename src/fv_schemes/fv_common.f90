!v Submodule file fv_common.smod
!
!  An implementation of the finite volume method

submodule(fv) fv_common
#include "ccs_macros.inc"
  use constants, only: add_mode, insert_mode, ndim
  use types, only: vector_values, matrix_values_spec, matrix_values, cell_locator, face_locator, &
                   neighbour_locator
  use vec, only: get_vector_data, restore_vector_data, create_vector_values

  use mat, only: create_matrix_values, set_matrix_values_spec_nrows, set_matrix_values_spec_ncols
  use utils, only: clear_entries, set_entry, set_row, set_col, set_values, set_mode, update
  use utils, only: debug_print, exit_print, str
  use meshing, only: count_neighbours, get_boundary_status, set_neighbour_location, &
                     get_local_index, get_global_index, get_volume, get_distance, &
                     set_face_location, get_face_area, get_face_normal, set_cell_location
  use boundary_conditions, only: get_bc_index
  use bc_constants

  implicit none

contains

  !> Computes fluxes and assign to matrix and RHS
  module subroutine compute_fluxes(phi, mf, mesh, component, cps, M, vec)
    class(field), intent(in) :: phi             
    class(field), intent(in) :: mf              
    type(ccs_mesh), intent(in) :: mesh          
    integer(ccs_int), intent(in) :: component   
    integer(ccs_int), intent(in) :: cps         
    class(ccs_matrix), intent(inout) :: M       
    class(ccs_vector), intent(inout) :: vec     

    integer(ccs_int) :: n_int_cells
    real(ccs_real), dimension(:), pointer :: mf_data

    associate (mf_values => mf%values)
      call dprint("CF: get mf")
      call get_vector_data(mf_values, mf_data)

      ! Loop over cells computing advection and diffusion fluxes
      n_int_cells = calc_matrix_nnz()
      call dprint("CF: compute coeffs")
      call compute_coeffs(phi, mf_data, mesh, component, n_int_cells, cps, M, vec)

      call dprint("CF: restore mf")
      call restore_vector_data(mf_values, mf_data)
    end associate

  end subroutine compute_fluxes

  !v Returns the number of entries per row that are non-zero
  !
  !  @note This assumes a square 2d grid
  pure function calc_matrix_nnz() result(nnz)
    integer(ccs_int) :: nnz   !< number of non-zero entries per row

    nnz = 5_ccs_int
  end function calc_matrix_nnz

  !> Computes the matrix coefficient for cells in the interior of the mesh
  subroutine compute_coeffs(phi, mf, mesh, component, n_int_cells, cps, M, b)
    class(field), intent(in) :: phi                !< scalar field structure
    real(ccs_real), dimension(:), intent(in) :: mf !< mass flux array defined at faces
    type(ccs_mesh), intent(in) :: mesh             !< Mesh structure
    integer(ccs_int), intent(in) :: component      !< integer indicating direction of velocity field component
    integer(ccs_int), intent(in) :: n_int_cells    !< number of cells in the interior of the mesh
    integer(ccs_int), intent(in) :: cps            !< number of cells per side
    class(ccs_matrix), intent(inout) :: M          !< Matrix structure being assigned
    class(ccs_vector), intent(inout) :: b          !< vector structure being assigned

    type(matrix_values_spec) :: mat_val_spec
    type(matrix_values) :: mat_coeffs
    type(vector_values) :: b_coeffs
    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb
    type(face_locator) :: loc_f
    integer(ccs_int) :: global_index_p, global_index_nb, index_p, index_nb
    integer(ccs_int) :: j
    integer(ccs_int) :: row, col
    integer(ccs_int) :: nnb
    real(ccs_real) :: face_area
    real(ccs_real) :: diff_coeff, diff_coeff_total
    real(ccs_real) :: adv_coeff, adv_coeff_total
    real(ccs_real) :: bc_value
    real(ccs_real), dimension(ndim) :: face_normal
    logical :: is_boundary

    integer(ccs_int) :: index_f

    real(ccs_real) :: sgn ! Sign indicating face orientation

    call set_matrix_values_spec_nrows(1_ccs_int, mat_val_spec)
    call set_matrix_values_spec_ncols(n_int_cells, mat_val_spec)
    call create_matrix_values(mat_val_spec, mat_coeffs)
    call set_mode(add_mode, mat_coeffs)

    call create_vector_values(n_int_cells, b_coeffs)
    call set_mode(add_mode, b_coeffs)

    do index_p = 1, mesh%nlocal
      call clear_entries(mat_coeffs)
      call clear_entries(b_coeffs)

      ! Calculate contribution from neighbours
      call set_cell_location(mesh, index_p, loc_p)
      call get_global_index(loc_p, global_index_p)
      call count_neighbours(loc_p, nnb)

      call set_row(global_index_p, mat_coeffs)

      adv_coeff_total = 0.0_ccs_real
      diff_coeff_total = 0.0_ccs_real
      do j = 1, nnb
        call set_neighbour_location(loc_p, j, loc_nb)
        call get_boundary_status(loc_nb, is_boundary)
        call set_face_location(mesh, index_p, j, loc_f)
        call get_face_normal(loc_f, face_normal)

        if (.not. is_boundary) then
          diff_coeff = calc_diffusion_coeff(index_p, j, mesh)

          call get_global_index(loc_nb, global_index_nb)
          call get_local_index(loc_nb, index_nb)

          call get_face_area(loc_f, face_area)
          call get_local_index(loc_f, index_f)

          ! XXX: Why won't Fortran interfaces distinguish on extended types...
          ! TODO: This will be expensive (in a tight loop) - investigate moving to a type-bound
          !       procedure (should also eliminate the type check).
          if (index_nb < index_p) then
            sgn = -1.0_ccs_real
          else
            sgn = 1.0_ccs_real
          end if
          select type (phi)
          type is (central_field)
            call calc_advection_coeff(phi, sgn * mf(index_f), 0, adv_coeff)
          type is (upwind_field)
            call calc_advection_coeff(phi, sgn * mf(index_f), 0, adv_coeff)
          class default
            call error_abort("Invalid velocity field discretisation.")
          end select

          ! XXX: we are relying on div(u)=0 => a_P = -sum_nb a_nb
          adv_coeff = adv_coeff * (sgn * mf(index_f) * face_area)

          call set_col(global_index_nb, mat_coeffs)
          call set_entry(adv_coeff + diff_coeff, mat_coeffs)
          adv_coeff_total = adv_coeff_total + adv_coeff
          diff_coeff_total = diff_coeff_total + diff_coeff
        else
          call get_local_index(loc_nb, index_nb)

          !call set_face_location(mesh, index_p, j, loc_f)
          call get_face_area(loc_f, face_area)
          call get_local_index(loc_f, index_f)

          diff_coeff = calc_diffusion_coeff(index_p, j, mesh)
          select type (phi)
          type is (central_field)
            call calc_advection_coeff(phi, mf(index_f), index_nb, adv_coeff)
          type is (upwind_field)
            call calc_advection_coeff(phi, mf(index_f), index_nb, adv_coeff)
          class default
            call error_abort("Invalid velocity field discretisation.")
          end select
          adv_coeff = adv_coeff * (mf(index_f) * face_area)

          !call calc_cell_coords(global_index_p, cps, row, col) ! XXX: don't think we need this anymore
          call compute_boundary_values(phi, component, index_nb, index_p, loc_p, loc_f, face_normal, bc_value)

          call set_row(global_index_p, b_coeffs)
          call set_entry(-(adv_coeff + diff_coeff) * bc_value, b_coeffs)

          call set_row(global_index_p, mat_coeffs)
          call set_col(global_index_p, mat_coeffs)
          call set_entry(-(adv_coeff + diff_coeff), mat_coeffs)
        end if
      end do

      call set_values(b_coeffs, b)

      call set_col(global_index_p, mat_coeffs)
      call set_entry(-(adv_coeff_total + diff_coeff_total), mat_coeffs)
      call set_values(mat_coeffs, M)
    end do

    deallocate (mat_coeffs%global_row_indices)
    deallocate (mat_coeffs%global_col_indices)
    deallocate (mat_coeffs%values)
  end subroutine compute_coeffs

  !> Computes the value of the scalar field on the boundary 
  subroutine compute_boundary_values(phi, component, index_nb, index_p, loc_p, loc_f, normal, bc_value, x_gradients, y_gradients, z_gradients)
    class(field), intent(in) :: phi                         !< the field for which boundary values are being computed
    integer(ccs_int), intent(in) :: component               !< integer indicating direction of velocity field component
    integer, intent(in) :: index_nb                         !< index of neighbour 
    integer, intent(in) :: index_p                          !< index of cell 
    type(cell_locator), intent(in) :: loc_p                 !< location of cell
    type(face_locator), intent(in) :: loc_f                 !< location of face
    real(ccs_real), dimension(ndim), intent(in) :: normal   !< boundary face normal direction
    real(ccs_real), intent(out) :: bc_value                 !< the value of the scalar field at the specified boundary
    real(ccs_real), dimension(:), optional, intent(in) :: x_gradients, y_gradients, z_gradients

    ! local variables
    integer(ccs_int) :: index_bc
    integer(ccs_int) :: i
    real(ccs_real), dimension(ndim) :: dx
    real(ccs_real), dimension(ndim) :: parallel_component_map
    real(ccs_real), dimension(ndim) :: phi_face_parallel_component
    real(ccs_real) :: phi_face_parallel_component_norm
    real(ccs_real) :: phi_face_parallel_component_portion
    real(ccs_real) :: normal_norm
    real(ccs_real), dimension(:), pointer :: phi_values

    call get_bc_index(phi, index_nb, index_bc)

    select case (phi%bcs%bc_types(index_bc))
    case (bc_type_dirichlet)
      bc_value = phi%bcs%values(index_bc)
    case (bc_type_extrapolate)
      call get_vector_data(phi%values, phi_values)
      call get_distance(loc_p, loc_f, dx)

      bc_value = phi_values(index_p) + (x_gradients(index_p) * dx(1) + y_gradients(index_p) * dx(2) + z_gradients(index_p) * dx(3))

      call restore_vector_data(phi%values, phi_values)
    case (bc_type_sym)
      select case(component)
      case (1)
        parallel_component_map = (/ 0, 1, 1 /)
      case (2)
        parallel_component_map = (/ 1, 0, 1 /)
      case (3)
        parallel_component_map = (/ 1, 1, 0 /)
      case default
        call error_abort("invalid component provided " // str(component))
      end select
      ! Only keep the components of phi that are parallel to the surface
      phi_face_parallel_component_norm = 0
      normal_norm = 0
      do i = 1, ndim
        phi_face_parallel_component(i) = parallel_component_map(i) * normal(i)
        phi_face_parallel_component_norm = phi_face_parallel_component_norm + phi_face_parallel_component(i) * phi_face_parallel_component(i)
        normal_norm = normal_norm + normal(i) * normal(i)
      end do
      phi_face_parallel_component_portion = sqrt(phi_face_parallel_component_norm/normal_norm)
      
      ! Get value of phi at boundary cell
      call get_vector_data(phi%values, phi_values)
      bc_value = phi_face_parallel_component_portion*phi_values(index_p)
    case default
      bc_value = 0.0_ccs_real
      call error_abort("unknown bc type " // str(phi%bcs%bc_types(index_bc)))
    end select
  end subroutine compute_boundary_values

  !> Sets the diffusion coefficient
  module function calc_diffusion_coeff(index_p, index_nb, mesh) result(coeff)
    integer(ccs_int), intent(in) :: index_p  !< the local cell index
    integer(ccs_int), intent(in) :: index_nb !< the local neigbouring cell index
    type(ccs_mesh), intent(in) :: mesh       !< the mesh structure
    real(ccs_real) :: coeff                  !< the diffusion coefficient

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

  !> Calculates mass flux across given face. Note: assumes rho = 1 and uniform grid
  module function calc_mass_flux(u, v, p, dpdx, dpdy, invAu, invAv, loc_f) result(flux)
    real(ccs_real), dimension(:), intent(in) :: u, v         !< arrays containing x, y velocities
    real(ccs_real), dimension(:), intent(in) :: p            !< array containing pressure
    real(ccs_real), dimension(:), intent(in) :: dpdx, dpdy   !< arrays containing pressure gradient in x and y
    real(ccs_real), dimension(:), intent(in) :: invAu, invAv !< arrays containing inverse momentum diagonal in x and y
    type(face_locator), intent(in) :: loc_f                  !< face locator

    real(ccs_real) :: flux                                   !< The flux across the boundary

    ! Local variables
    logical :: is_boundary                         ! Boundary indicator
    type(cell_locator) :: loc_p                    ! Primary cell locator
    type(neighbour_locator) :: loc_nb              ! Neighbour cell locator
    integer(ccs_int) :: index_nb                   ! Neighbour cell index
    real(ccs_real) :: flux_corr                    ! Flux correction
    real(ccs_real), dimension(ndim) :: dx          ! Cell-cell distance
    real(ccs_real) :: dxmag                        ! Cell-cell distance magnitude
    real(ccs_real), dimension(ndim) :: face_normal ! (local) face-normal array
    real(ccs_real) :: Vp                           ! Primary cell volume
    real(ccs_real) :: V_nb                         ! Neighbour cell volume
    real(ccs_real) :: Vf                           ! Face "volume"
    real(ccs_real) :: invAp                        ! Primary cell inverse momentum coefficient
    real(ccs_real) :: invA_nb                      ! Neighbour cell inverse momentum coefficient
    real(ccs_real) :: invAf                        ! Face inverse momentum coefficient

    call get_boundary_status(loc_f, is_boundary)
    if (.not. is_boundary) then
      associate (mesh => loc_f%mesh, &
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
          ! XXX: making convention to point from low to high cell.
          flux = -flux
        end if
      end associate
    else
      ! TODO: Write more general implementation handling BCs
      flux = 0.0_ccs_real ! XXX: hardcoded zero-flux BC
    end if

  end function calc_mass_flux

  !> Calculates the row and column indices from flattened vector index. Assumes square mesh
  module subroutine calc_cell_coords(index, cps, row, col)
    integer(ccs_int), intent(in) :: index !< cell index
    integer(ccs_int), intent(in) :: cps   !< number of cells per side
    integer(ccs_int), intent(out) :: row  !< cell row within mesh
    integer(ccs_int), intent(out) :: col  !< cell column within mesh

    col = modulo(index - 1, cps) + 1
    row = (index - 1) / cps + 1
  end subroutine calc_cell_coords

  !v Performs an update of the gradients of a field.
  !  @note This will perform a parallel update of the gradient fields to ensure halo cells are
  !  correctly updated on other PEs. @endnote
  module subroutine update_gradient(mesh, phi)

    type(ccs_mesh), intent(in) :: mesh !< the mesh
    class(field), intent(inout) :: phi !< the field whose gradients we want to update

    real(ccs_real), dimension(:), pointer :: x_gradients_data, y_gradients_data, z_gradients_data
    real(ccs_real), dimension(:), allocatable :: x_gradients_old, y_gradients_old, z_gradients_old

    integer(ccs_real) :: i

    call get_vector_data(phi%x_gradients, x_gradients_data)
    call get_vector_data(phi%y_gradients, y_gradients_data)
    call get_vector_data(phi%z_gradients, z_gradients_data)

    associate (ntotal => mesh%ntotal)
      allocate (x_gradients_old(ntotal))
      allocate (y_gradients_old(ntotal))
      allocate (z_gradients_old(ntotal))
      do i = 1, ntotal
        x_gradients_old(i) = x_gradients_data(i)
        y_gradients_old(i) = y_gradients_data(i)
        z_gradients_old(i) = z_gradients_data(i)
      end do
    end associate

    call restore_vector_data(phi%x_gradients, x_gradients_data)
    call restore_vector_data(phi%y_gradients, y_gradients_data)
    call restore_vector_data(phi%z_gradients, z_gradients_data)

    call dprint("Compute x gradient")
    call update_gradient_component(mesh, 1, phi, x_gradients_old, y_gradients_old, z_gradients_old, phi%x_gradients)
    call update(phi%x_gradients) ! XXX: opportunity to overlap update with later compute (begin/compute/end)
    call dprint("Compute y gradient")
    call update_gradient_component(mesh, 2, phi, x_gradients_old, y_gradients_old, z_gradients_old, phi%y_gradients)
    call update(phi%y_gradients) ! XXX: opportunity to overlap update with later compute (begin/compute/end)
    call dprint("Compute z gradient")
    call update_gradient_component(mesh, 3, phi, x_gradients_old, y_gradients_old, z_gradients_old, phi%z_gradients)
    call update(phi%z_gradients) ! XXX: opportunity to overlap update with later compute (begin/compute/end)

    deallocate (x_gradients_old)
    deallocate (y_gradients_old)
    deallocate (z_gradients_old)

  end subroutine update_gradient

  !> Helper subroutine to calculate a gradient component at a time.
  subroutine update_gradient_component(mesh, component, phi, x_gradients_old, y_gradients_old, z_gradients_old, gradients)

    type(ccs_mesh), intent(in) :: mesh !< the mesh
    integer(ccs_int), intent(in) :: component !< which vector component (i.e. direction) to update?
    class(field), intent(in) :: phi !< the field whose gradient we want to compute
    real(ccs_real), dimension(:), intent(in) :: x_gradients_old
    real(ccs_real), dimension(:), intent(in) :: y_gradients_old
    real(ccs_real), dimension(:), intent(in) :: z_gradients_old
    class(ccs_vector), intent(inout) :: gradients !< a cell-centred array of the gradient

    type(vector_values) :: grad_values
    real(ccs_real), dimension(:), pointer :: phi_data
    real(ccs_real) :: grad

    integer(ccs_int) :: index_p
    integer(ccs_int) :: j
    type(cell_locator) :: loc_p
    type(face_locator) :: loc_f
    type(neighbour_locator) :: loc_nb

    integer(ccs_int) :: nnb
    integer(ccs_int) :: index_nb

    real(ccs_real) :: phif

    logical :: is_boundary

    real(ccs_real) :: face_area
    real(ccs_real), dimension(ndim) :: face_norm

    real(ccs_real) :: V
    integer(ccs_int) :: global_index_p

    real(ccs_real), dimension(ndim) :: dx

    call create_vector_values(1_ccs_int, grad_values)
    call set_mode(insert_mode, grad_values)

    call get_vector_data(phi%values, phi_data)

    do index_p = 1, mesh%nlocal
      call clear_entries(grad_values)

      grad = 0.0_ccs_int

      call set_cell_location(mesh, index_p, loc_p)
      call count_neighbours(loc_p, nnb)
      do j = 1, nnb
        call set_face_location(mesh, index_p, j, loc_f)
        call get_boundary_status(loc_f, is_boundary)
        call get_face_area(loc_f, face_area)
        call get_face_normal(loc_f, face_norm)

        call set_neighbour_location(loc_p, j, loc_nb)
        call get_local_index(loc_nb, index_nb)
        if (.not. is_boundary) then
          phif = 0.5_ccs_real * (phi_data(index_p) + phi_data(index_nb)) ! XXX: Need to do proper interpolation
        else
          ! XXX: Add boundary condition treatment
          !call get_distance(loc_p, loc_f, dx)
          !phif = phi_data(index_p) + (x_gradients_old(index_p) * dx(1) + y_gradients_old(index_p) * dx(2) + z_gradients_old(index_p) * dx(3))
          call compute_boundary_values(phi, component, index_nb, index_p, loc_p, loc_f, face_norm, phif, x_gradients_old, y_gradients_old, z_gradients_old)
        end if

        grad = grad + phif * (face_area * face_norm(component))
      end do

      call get_volume(loc_p, V)
      grad = grad / V

      call get_global_index(loc_p, global_index_p)
      call set_row(global_index_p, grad_values)
      call set_entry(grad, grad_values)
      call set_values(grad_values, gradients)
    end do

    call restore_vector_data(phi%values, phi_data)

    deallocate (grad_values%global_indices)
    deallocate (grad_values%values)

  end subroutine update_gradient_component

end submodule fv_common
