!v Submodule file fv_common.smod
!
!  An implementation of the finite volume method
submodule(fv) fv_common
#include "ccs_macros.inc"
  use constants, only: insert_mode, add_mode
  use types, only: vector_values, matrix_values_spec, matrix_values, neighbour_locator, bc_profile, field, field_ptr
  use vec, only: get_vector_data, restore_vector_data, &
                 get_vector_data_readonly, restore_vector_data_readonly, &
                 create_vector_values, begin_ghost_update_vector, end_ghost_update_vector

  use mat, only: create_matrix_values, set_matrix_values_spec_nrows, set_matrix_values_spec_ncols
  use utils, only: clear_entries, set_entry, set_row, set_col, set_values, set_mode, update
  use utils, only: debug_print, exit_print, str
  use fv_equations, only: momentum_equation
  use fv_kernels, only: advection_kernel, create_advection_kernel
  use meshing, only: count_neighbours, get_boundary_status, create_neighbour_locator, &
                     get_local_index, get_global_index, get_volume, get_distance, &
                     create_face_locator, get_face_area, get_face_normal, create_cell_locator, &
                     get_local_num_cells, get_face_interpolation, &
                     get_max_faces, get_centre
  use boundary_conditions, only: get_bc_index
  use profiler, only: profiler_begin_region, profiler_end_region
  use bc_constants
  use error_codes

  implicit none

contains

  !> Computes fluxes and assign to matrix and RHS.
  !> Uses the mesh-wide maximum neighbour face count to size owner+neighbour stencil buffers.
  module subroutine compute_fluxes(phi, mf, viscosity, density, component, M, vec)
    class(field), intent(inout) :: phi
    class(field), intent(inout) :: mf
    class(field), intent(inout) :: viscosity
    class(field), intent(inout) :: density
    integer(ccs_int), intent(in) :: component
    class(ccs_matrix), intent(inout) :: M
    class(ccs_vector), intent(inout) :: vec

    integer(ccs_int) :: max_nb_faces !< Maximum neighbour faces per cell (mesh-wide upper bound)
    integer(ccs_int) :: n_stencil_cells !< Owner cell plus its maximum neighbour faces
    real(ccs_real), dimension(:), pointer :: mf_data, viscosity_data, density_data

    associate (mf_values => mf%values)
      call dprint("CF: get mf")
      call get_vector_data(mf_values, mf_data)
      call get_vector_data(viscosity%values, viscosity_data)
      call get_vector_data(density%values, density_data)

      ! Loop over cells computing advection and diffusion fluxes
      call get_max_faces(max_nb_faces)
      n_stencil_cells = max_nb_faces + 1 ! 1 neighbour per face + central cell
      call dprint("CF: compute coeffs")
      call compute_coeffs(phi, mf_data, viscosity_data, density_data, component, n_stencil_cells, M, vec)

      call dprint("CF: restore mf")
      call restore_vector_data(mf_values, mf_data)
      call restore_vector_data(viscosity%values, viscosity_data)
      call restore_vector_data(density%values, density_data)
    end associate

  end subroutine compute_fluxes

  !> Computes the matrix coefficient for cells in the interior of the mesh
  subroutine compute_coeffs(phi, mf, visc, dens, component, n_stencil_cells, M, b)

    class(field), intent(inout) :: phi !< scalar field structure
    real(ccs_real), dimension(:), target, intent(in) :: mf   !< mass flux array defined at faces
    real(ccs_real), dimension(:), target, intent(in) :: visc !< viscosity
    real(ccs_real), dimension(:), target, intent(in) :: dens !< density
    integer(ccs_int), intent(in) :: component                !< integer indicating direction of velocity field component
    integer(ccs_int), intent(in) :: n_stencil_cells          !< Owner cell plus the maximum neighbour faces
    class(ccs_matrix), intent(inout) :: M                    !< equation system matrix
    class(ccs_vector), intent(inout) :: b                    !< RHS vector

    type(matrix_values_spec) :: mat_val_spec
    type(matrix_values) :: mat_coeffs
    type(vector_values) :: b_coeffs
    type(cell_locator) :: loc_p
    type(momentum_equation) :: mom_eq
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: index_p
    integer(ccs_int) :: max_nb_faces
    integer(ccs_int) :: max_kernel_width
    integer(ccs_int) :: n_cells_eff
    class(advection_kernel), allocatable :: kernel

    max_nb_faces = max(0_ccs_int, n_stencil_cells - 1_ccs_int)

    kernel = create_advection_kernel(phi)
    max_kernel_width = max(mom_eq%diff_kernel%get_width(), kernel%get_width())
    n_cells_eff = 1_ccs_int + max_nb_faces * max_kernel_width

    call set_matrix_values_spec_nrows(1_ccs_int, mat_val_spec)
    call set_matrix_values_spec_ncols(n_cells_eff, mat_val_spec)
    call create_matrix_values(mat_val_spec, mat_coeffs)
    call set_mode(add_mode, mat_coeffs)

    call create_vector_values(n_cells_eff, b_coeffs)
    call set_mode(add_mode, b_coeffs)

    call mom_eq%set_advection(kernel)
    call mom_eq%init(max_nb_faces, mf, visc, dens, component)

    call get_local_num_cells(local_num_cells)
    do index_p = 1, local_num_cells
      call clear_entries(mat_coeffs)
      call clear_entries(b_coeffs)
      call create_cell_locator(index_p, loc_p)

      call mom_eq%gather(phi, loc_p)
      call mom_eq%apply(mat_coeffs, b_coeffs)

      call set_values(b_coeffs, b)
      call set_values(mat_coeffs, M)
    end do

    call update(M)
    call update(b)

    deallocate(mat_coeffs%global_row_indices)
    deallocate(mat_coeffs%global_col_indices)
    deallocate(mat_coeffs%values)
  end subroutine compute_coeffs

  !> Computes the value of the scalar field on the boundary
  pure module subroutine compute_boundary_values(phi, component, loc_p, loc_f, normal, bc_value)

    class(field), intent(in) :: phi                       !< the field for which boundary values are being computed
    integer(ccs_int), intent(in) :: component             !< integer indicating direction of velocity field component
    type(cell_locator), intent(in) :: loc_p               !< location of cell
    type(face_locator), intent(in) :: loc_f               !< location of face
    real(ccs_real), dimension(ndim), intent(in) :: normal !< boundary face normal direction
    real(ccs_real), intent(out) :: bc_value               !< the boundary value

    real(ccs_real) :: a !< The diagonal coeff (implicit component)
    real(ccs_real) :: b !< The RHS value (explicit component)
    integer(ccs_int) :: index_p

    call compute_boundary_coeffs(phi, component, loc_p, loc_f, normal, a, b)

    call get_local_index(loc_p, index_p)
    bc_value = 0.5_ccs_real * (phi%values_ro(index_p) + (b + a * phi%values_ro(index_p)))

  end subroutine compute_boundary_values

  !> Compute the coefficients of the boundary condition
  pure module subroutine compute_boundary_coeffs(phi, component, loc_p, loc_f, normal, a, b)

    class(field), intent(in) :: phi                       !< the field for which boundary values are being computed
    integer(ccs_int), intent(in) :: component             !< integer indicating direction of velocity field component
    type(cell_locator), intent(in) :: loc_p               !< location of cell
    type(face_locator), intent(in) :: loc_f               !< location of face
    real(ccs_real), dimension(ndim), intent(in) :: normal !< boundary face normal direction
    real(ccs_real), intent(out) :: a                      !< The diagonal coeff (implicit)
    real(ccs_real), intent(out) :: b                      !< The RHS entry (explicit)

    ! local variables
    integer(ccs_int) :: index_bc
    integer(ccs_int) :: index_nb
    integer(ccs_int) :: index_p
    type(neighbour_locator) :: loc_nb
    ! integer(ccs_int) :: i
    real(ccs_real), dimension(ndim) :: dx
    real(ccs_real), dimension(ndim) :: x
    ! real(ccs_real), dimension(ndim) :: parallel_component_map
    ! real(ccs_real), dimension(ndim) :: phi_face_parallel_component
    ! real(ccs_real) :: phi_face_parallel_component_norm
    ! real(ccs_real) :: phi_face_parallel_component_portion
    ! real(ccs_real) :: normal_norm
    real(ccs_real) :: dxmag
    real(ccs_real) :: bc_value

    associate(foo=> component, bar=>normal)
      end associate
    call get_local_index(loc_p, index_p)
    call create_neighbour_locator(loc_p, loc_f%cell_face_ctr, loc_nb)
    call get_local_index(loc_nb, index_nb)
    call get_bc_index(phi, index_nb, index_bc)

    select case (phi%bcs%bc_types(index_bc))
    case (bc_type_dirichlet)
      a = -1.0_ccs_real
      b = 2.0_ccs_real * phi%bcs%values(index_bc)
    case (bc_type_extrapolate)
      call get_distance(loc_p, loc_f, dx)

      a = 1.0_ccs_real
b = 2.0_ccs_real * (phi%x_gradients_ro(index_p) * dx(1) + phi%y_gradients_ro(index_p) * dx(2) + phi%z_gradients_ro(index_p) * dx(3))
    ! case (bc_type_sym)  ! XXX: Make sure this works as intended for symmetric BC.
    !   select case (component)
    !   case (0)
    !     parallel_component_map = [1, 1, 1]
    !   case (1)
    !     parallel_component_map = [0, 1, 1]
    !   case (2)
    !     parallel_component_map = [1, 0, 1]
    !   case (3)
    !     parallel_component_map = [1, 1, 0]
    !   case default
    !     error stop invalid_component ! Invalid component provided
    !   end select
    !   ! Only keep the components of phi that are parallel to the surface
    !   phi_face_parallel_component_norm = 0
    !   normal_norm = 0
    !   do i = 1, ndim
    !     phi_face_parallel_component(i) = parallel_component_map(i) * normal(i)
    !     phi_face_parallel_component_norm = phi_face_parallel_component_norm + &
    !                                        phi_face_parallel_component(i) * phi_face_parallel_component(i)
    !     normal_norm = normal_norm + normal(i) * normal(i)
    !   end do
    !   phi_face_parallel_component_portion = sqrt(phi_face_parallel_component_norm / normal_norm)

    !   ! Get value of phi at boundary cell
    !   a = phi_face_parallel_component_portion
    !   b = 0.0_ccs_real
    case (bc_type_neumann)
      call get_distance(loc_p, loc_f, dx)
      dxmag = norm2(dx)

      a = 1.0_ccs_real
      b = (2.0_ccs_real * dxmag) * phi%bcs%values(index_bc)
    case (bc_type_profile)
      call get_centre(loc_f, x)
      if (allocated(phi%bcs%profiles(index_bc)%centre)) then
        call get_value_from_bc_profile(x, phi%bcs%profiles(index_bc), bc_value)
      else
        bc_value = 0.0_ccs_real
      end if

      a = -1.0_ccs_real
      b = 2.0_ccs_real * bc_value
    case default
      ! Set coefficients to cause divergence
      ! Prevents "unused variable" compiler errors
      a = 0.0_ccs_real
      b = huge(1.0_ccs_real)

      error stop unknown_bc_type ! Unknown BC type
    end select

  end subroutine

  !> Linear interpolate of BC profile
  pure module subroutine get_value_from_bc_profile(x, profile, bc_value)
    real(ccs_real), dimension(:), intent(in) :: x
    type(bc_profile), intent(in) :: profile
    real(ccs_real), intent(out) :: bc_value
    integer(ccs_int) :: n, i
    real(ccs_real) :: r
    real(ccs_real) :: coeff

    r = norm2(x(:) - profile%centre(:))

    n = size(profile%coordinates)

    bc_value = profile%values(n)
    if (r .le. profile%coordinates(1)) then
      bc_value = profile%values(1)
      return
    end if

    do i = 1, n - 1
      if (r .lt. profile%coordinates(i + 1)) then
        coeff = (r - profile%coordinates(i)) / (profile%coordinates(i + 1) - profile%coordinates(i))
        bc_value = (1 - coeff) * profile%values(i) + coeff * profile%values(i + 1)
        return
      end if
    end do

  end subroutine

  !> Sets the diffusion coefficient
  pure module subroutine calc_diffusion_coeff(index_p, index_nb, enable_cell_corrections, visc_p, visc_nb, dens_p, dens_nb, SchmidtNo, coeff) 
    integer(ccs_int), intent(in) :: index_p  !< the local cell index
    integer(ccs_int), intent(in) :: index_nb !< the local neigbouring cell index
    logical, intent(in) :: enable_cell_corrections !< Whether or not cell shape corrections are used
    real(ccs_real), intent(out) :: coeff                  !< the diffusion coefficient
    real(ccs_real), intent(in) :: visc_p, visc_nb !< viscosity
    real(ccs_real), intent(in) :: SchmidtNo
    real(ccs_real), intent(in) :: dens_p, dens_nb !< density
    type(face_locator) :: loc_f
    real(ccs_real) :: face_area
    real(ccs_real) :: diffusion_factor
    logical :: is_boundary
    real(ccs_real), dimension(ndim) :: dx
    real(ccs_real), dimension(ndim) :: n
    real(ccs_real), dimension(ndim) :: x_p
    real(ccs_real), dimension(ndim) :: x_nb
    real(ccs_real), dimension(ndim) :: x_f
    real(ccs_real) :: dxmag
    real(ccs_real) :: dx_orth
    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb
    real(ccs_real) :: visc_avg !< average viscosity
    real(ccs_real) :: dens_avg !< average density
    real(ccs_real) :: interpolation_factor

    call create_face_locator(index_p, index_nb, loc_f)
    call get_face_area(loc_f, face_area)
    call get_boundary_status(loc_f, is_boundary)

    call create_cell_locator(index_p, loc_p)
    if (.not. is_boundary) then
      call get_face_interpolation(loc_f, interpolation_factor)
      call create_neighbour_locator(loc_p, index_nb, loc_nb)

      if (enable_cell_corrections) then
        call get_face_normal(loc_f, n)

        call get_centre(loc_p, x_p)
        call get_centre(loc_nb, x_nb)
        call get_centre(loc_f, x_f)

        ! see math below, but it works because ||n||=1 and points outwards
        !rnb_k_prime = x_f + a*n
        !rp_prime = x_f - a*n
        !dx = rnb_k_prime - rp_prime
        !dxmag = norm2(dx)
        dx_orth = min(dot_product(x_f - x_p, n), dot_product(x_nb - x_f, n))
        dxmag = 2.0_ccs_real * dx_orth
      else
        call get_distance(loc_p, loc_nb, dx)
        dxmag = norm2(dx)
      end if
    else
      call get_distance(loc_p, loc_f, dx)
      dxmag = norm2(dx)
    end if

    if (.not. is_boundary) then
      visc_avg = (interpolation_factor * visc_p) + ((1.0_ccs_real - interpolation_factor) * visc_nb)
      dens_avg = (interpolation_factor * dens_p) + ((1.0_ccs_real - interpolation_factor) * dens_nb)
      diffusion_factor = visc_avg / (dens_avg * SchmidtNo)
    else
      diffusion_factor = visc_p / (dens_p * SchmidtNo)
    end if
    
    coeff = face_area * diffusion_factor / dxmag
  end subroutine calc_diffusion_coeff

  !> Interpolate field to face center from cell center, applied gradient correction (if enabled in the field
  ! spec) using Ferziger & Peric 4th ed, sec 9.7.1
  pure subroutine interpolate_field_to_face(phi, loc_f, face_value, face_correction_only)

    class(field), intent(inout) :: phi
    type(face_locator), intent(in) :: loc_f                         !< face locator
    real(ccs_real), intent(out) :: face_value
    real(ccs_real), optional, intent(out) :: face_correction_only

    real(ccs_real) :: face_correction
    type(cell_locator) :: loc_p                    ! Primary cell locator
    type(neighbour_locator) :: loc_nb              ! Neighbour cell locator
    integer(ccs_int) :: index_nb                   ! Neighbour cell index
    real(ccs_real), dimension(ndim) :: n           ! (local) face-normal array
    real(ccs_real), dimension(ndim) :: grad_phi_p, grad_phi_nb
    real(ccs_real), dimension(ndim) :: x_nb, x_p, x_f, x_nb_prime, x_p_prime
    real(ccs_real) :: interpol_factor, dx_orth

    associate (index_p => loc_f%index_p, &
               j => loc_f%cell_face_ctr)

      call create_cell_locator(index_p, loc_p)
      call create_neighbour_locator(loc_p, j, loc_nb)
      call get_local_index(loc_nb, index_nb)

      if (phi%enable_cell_corrections) then
        call get_face_normal(loc_f, n)
        call get_centre(loc_p, x_p)
        call get_centre(loc_nb, x_nb)
        call get_centre(loc_f, x_f)

        dx_orth = min(dot_product(x_f - x_p, n), dot_product(x_nb - x_f, n))
        x_nb_prime = x_f + dx_orth * n
        x_p_prime = x_f - dx_orth * n

        grad_phi_p = [phi%x_gradients_ro(index_p), phi%y_gradients_ro(index_p), phi%z_gradients_ro(index_p)]
        grad_phi_nb = [phi%x_gradients_ro(index_nb), phi%y_gradients_ro(index_nb), phi%z_gradients_ro(index_nb)]

        face_correction = 0.5_ccs_real * (dot_product(grad_phi_p, x_p_prime - x_p) + dot_product(grad_phi_nb, x_nb_prime - x_nb))
        face_value = 0.5_ccs_real * (phi%values_ro(index_p) + phi%values_ro(index_nb)) + face_correction
      else
        call get_face_interpolation(loc_f, interpol_factor)
        face_correction = 0.0_ccs_real
        face_value = (interpol_factor * phi%values_ro(index_p) + (1.0_ccs_real - interpol_factor) * phi%values_ro(index_nb))
      end if

      if (present(face_correction_only)) then
        face_correction_only = face_correction
      end if

    end associate

  end subroutine

  !> Calculates mass flux across given face. Note: assumes rho = 1
  module function calc_mass_flux_uvw(u_field, v_field, w_field, p, dpdx, dpdy, dpdz, invA, &
                                     loc_f, enable_cell_corrections) result(flux)
    class(field), intent(inout) :: u_field
    class(field), intent(inout) :: v_field
    class(field), intent(inout) :: w_field
    real(ccs_real), dimension(:), intent(in) :: p                !< array containing pressure
    real(ccs_real), dimension(:), intent(in) :: dpdx, dpdy, dpdz !< arrays containing pressure gradient in x, y and z
    real(ccs_real), dimension(:), intent(in) :: invA             !< array containing inverse momentum diagonal
    type(face_locator), intent(in) :: loc_f                      !< face locator
    logical, intent(in) :: enable_cell_corrections

    real(ccs_real) :: flux                                       !< The flux across the boundary

    ! Local variables
    logical :: is_boundary                         ! Boundary indicator
    type(cell_locator) :: loc_p                    ! Primary cell locator
    type(neighbour_locator) :: loc_nb              ! Neighbour cell locator
    integer(ccs_int) :: index_nb                   ! Neighbour cell index
    real(ccs_real) :: flux_corr                    ! Flux correction
    real(ccs_real), dimension(ndim) :: n           ! (local) face-normal array
    real(ccs_real) :: u_bc, v_bc, w_bc             ! values of u, v and w at boundary
    integer(ccs_int), parameter :: x_direction = 1, y_direction = 2, z_direction = 3
    real(ccs_real) :: flux_x, flux_y, flux_z

    call get_boundary_status(loc_f, is_boundary)

    associate (index_p => loc_f%index_p, &
               j => loc_f%cell_face_ctr)

      call create_cell_locator(index_p, loc_p)
      call create_neighbour_locator(loc_p, j, loc_nb)
      call get_local_index(loc_nb, index_nb)

      call get_face_normal(loc_f, n)
      ! XXX: this is likely expensive inside a loop...
      if (.not. is_boundary) then

        call interpolate_field_to_face(u_field, loc_f, flux_x)
        call interpolate_field_to_face(v_field, loc_f, flux_y)
        call interpolate_field_to_face(w_field, loc_f, flux_z)
        flux = dot_product([flux_x, flux_y, flux_z], n)

        if (index_p > index_nb) then
          ! XXX: making convention to point from low to high cell.
          flux = -flux
        end if

        flux_corr = calc_mass_flux(p, dpdx, dpdy, dpdz, invA, loc_f, enable_cell_corrections)
        flux = flux + flux_corr
      else
        call compute_boundary_values(u_field, x_direction, loc_p, loc_f, n, u_bc)
        call compute_boundary_values(v_field, y_direction, loc_p, loc_f, n, v_bc)
        call compute_boundary_values(w_field, z_direction, loc_p, loc_f, n, w_bc)
        flux = dot_product([u_bc, v_bc, w_bc], n)
      end if
    end associate

  end function calc_mass_flux_uvw

  ! Computes Rhie-Chow correction
  pure module function calc_mass_flux_no_uvw(p, dpdx, dpdy, dpdz, invA, loc_f, enable_cell_corrections) result(flux)
    real(ccs_real), dimension(:), intent(in) :: p                !< array containing pressure
    real(ccs_real), dimension(:), intent(in) :: dpdx, dpdy, dpdz !< arrays containing pressure gradient in x, y and z
    real(ccs_real), dimension(:), intent(in) :: invA             !< array containing inverse momentum diagonal
    type(face_locator), intent(in) :: loc_f                      !< face locator
    logical, intent(in) :: enable_cell_corrections

    real(ccs_real) :: flux                         !< The flux across the boundary

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
    integer(ccs_int), parameter :: x_direction = 1, y_direction = 2, z_direction = 3

    real(ccs_real) :: interpol_factor
    real(ccs_real), dimension(ndim) :: grad_phi_p
    real(ccs_real), dimension(ndim) :: grad_phi_nb
    real(ccs_real), dimension(ndim) :: x_nb, x_p, x_f, rnb_k_prime, rp_prime
    real(ccs_real), dimension(ndim) :: n
    real(ccs_real) :: dx_orth

    call get_boundary_status(loc_f, is_boundary)

    associate (index_p => loc_f%index_p, &
               j => loc_f%cell_face_ctr)

      call create_cell_locator(index_p, loc_p)
      call create_neighbour_locator(loc_p, j, loc_nb)
      call get_local_index(loc_nb, index_nb)

      call get_face_normal(loc_f, face_normal)
      if (.not. is_boundary) then

        !
        ! Rhie-Chow correction from Ferziger & Peric
        !
        call get_face_interpolation(loc_f, interpol_factor)
        call get_face_normal(loc_f, face_normal)

        grad_phi_p = [dpdx(index_p), dpdy(index_p), dpdz(index_p)]
        grad_phi_nb = [dpdx(index_nb), dpdy(index_nb), dpdz(index_nb)]

        if (enable_cell_corrections) then
          ! Cell excentricity/non-orthogonality corrections (Ferziger & Peric 4th ed, sec 9.8, p317, eq9.67 and 9.66)
          call get_face_normal(loc_f, n)
          call get_centre(loc_p, x_p)
          call get_centre(loc_nb, x_nb)
          call get_centre(loc_f, x_f)

          dx_orth = min(dot_product(x_f - x_p, n), dot_product(x_nb - x_f, n))
          rnb_k_prime = x_f + dx_orth * n
          rp_prime = x_f - dx_orth * n
          !dx = rnb_k_prime - rp_prime
          !dxmag = norm2(dx)
          dxmag = 2.0_ccs_real * dx_orth

          ! eq 9.66
          flux_corr = (p(index_nb) - p(index_p)) + (dot_product(grad_phi_nb, rnb_k_prime - x_nb) - dot_product(grad_phi_p, rp_prime - x_p))

          ! interpolated pressure gradient at the face (i.e.)
          flux_corr = flux_corr - 0.5_ccs_real * dot_product((grad_phi_p + grad_phi_nb), x_nb - x_p)

          flux_corr = -flux_corr / dxmag
        else
          call get_distance(loc_p, loc_nb, dx)
          dxmag = norm2(dx)
          flux_corr = -(p(index_nb) - p(index_p)) / dxmag
          flux_corr = flux_corr + dot_product(interpol_factor * grad_phi_p + (1.0_ccs_real - interpol_factor) * grad_phi_nb, face_normal)
        end if

        call get_volume(loc_p, Vp)
        call get_volume(loc_nb, V_nb)
        Vf = interpol_factor * Vp + (1.0_ccs_real - interpol_factor) * V_nb

        ! This is probably not quite right ...
        invAp = invA(index_p)
        invA_nb = invA(index_nb)
        invAf = interpol_factor * invAp + (1.0_ccs_real - interpol_factor) * invA_nb

        flux_corr = (Vf * invAf) * flux_corr

        if (index_p > index_nb) then
          ! XXX: making convention to point from low to high cell.
          flux_corr = -flux_corr
        end if

        ! Apply correction
        flux = flux_corr
      else
        flux = 0.0_ccs_real
      end if
    end associate
  end function calc_mass_flux_no_uvw

  !> Calculates the row and column indices from flattened vector index. Assumes square mesh
  pure module subroutine calc_cell_coords(index, cps, row, col)
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
  module subroutine update_gradient_field(phi)

    class(field), target, intent(inout) :: phi !< the field whose gradients we want to update
    type(field_ptr), dimension(1) :: fields

    fields(1)%ptr => phi
    call update_gradient_fields(fields)

  end subroutine update_gradient_field

  !v Performs an update of the gradients for several fields, overlapping communication between fields.
  module subroutine update_gradient_fields(fields)

    type(field_ptr), dimension(:), intent(inout) :: fields !< the fields whose gradients we want to update
    real(ccs_real), dimension(:, :), allocatable :: gradients
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: i
    integer(ccs_int) :: nfields

    call profiler_begin_region("Compute gradient")

    nfields = size(fields)
    if (nfields == 0) then
      call profiler_end_region("Compute gradient")
      return
    end if

    call get_local_num_cells(local_num_cells)
    allocate (gradients(local_num_cells, ndim))

    do i = 1, nfields
      call compute_gradients(fields(i)%ptr, gradients)
      call set_gradients(fields(i)%ptr, gradients)
      call start_gradient_halo(fields(i)%ptr)
    end do

    do i = 1, nfields
      call finish_gradient_halo(fields(i)%ptr)
    end do

    call profiler_end_region("Compute gradient")

  end subroutine update_gradient_fields

  subroutine compute_gradients(phi, gradients)
    class(field), intent(in) :: phi
    real(ccs_real), dimension(:, :), intent(out) :: gradients

    integer(ccs_int) :: i
    integer(ccs_int) :: local_num_cells
    type(cell_locator) :: loc_p
    real(ccs_real), dimension(ndim) :: grad_p

    local_num_cells = size(gradients, 1)
    do i = 1, local_num_cells
      call create_cell_locator(i, loc_p)
      call compute_gradient_at_point(phi, loc_p, grad_p)
      gradients(i, :) = grad_p(:)
    end do
  end subroutine compute_gradients

  subroutine set_gradients(phi, gradients)
    class(field), intent(inout) :: phi
    real(ccs_real), dimension(:, :), intent(in) :: gradients

    real(ccs_real), dimension(:), pointer :: x_gradient_data
    real(ccs_real), dimension(:), pointer :: y_gradient_data
    real(ccs_real), dimension(:), pointer :: z_gradient_data
    integer(ccs_int) :: local_num_cells

    local_num_cells = size(gradients, 1)

    call get_vector_data(phi%x_gradients, x_gradient_data)
    call get_vector_data(phi%y_gradients, y_gradient_data)
    call get_vector_data(phi%z_gradients, z_gradient_data)

    x_gradient_data(1:local_num_cells) = gradients(:, 1)
    y_gradient_data(1:local_num_cells) = gradients(:, 2)
    z_gradient_data(1:local_num_cells) = gradients(:, 3)

    call restore_vector_data(phi%x_gradients, x_gradient_data)
    call restore_vector_data(phi%y_gradients, y_gradient_data)
    call restore_vector_data(phi%z_gradients, z_gradient_data)
  end subroutine set_gradients

  subroutine start_gradient_halo(phi)
    class(field), intent(inout) :: phi

    call begin_ghost_update_vector(phi%x_gradients)
    call begin_ghost_update_vector(phi%y_gradients)
    call begin_ghost_update_vector(phi%z_gradients)
  end subroutine start_gradient_halo

  subroutine finish_gradient_halo(phi)
    class(field), intent(inout) :: phi

    call end_ghost_update_vector(phi%x_gradients)
    call end_ghost_update_vector(phi%y_gradients)
    call end_ghost_update_vector(phi%z_gradients)
  end subroutine finish_gradient_halo

  !> Helper subroutine to calculate a gradient at a point (cell centre)
  pure subroutine compute_gradient_at_point(phi, loc_p, gradients)

    class(field), intent(in) :: phi         !< the field whose gradient we want to compute
    type(cell_locator), intent(in) :: loc_p !< locator of the current cell gradient is being evaluated in
    real(ccs_real), dimension(3), intent(out) :: gradients !< a cell-centred array of the gradient

    integer(ccs_int) :: index_p
    integer(ccs_int) :: j
    type(face_locator) :: loc_f
    type(neighbour_locator) :: loc_nb

    integer(ccs_int) :: nnb
    integer(ccs_int) :: index_nb

    real(ccs_real) :: phif
    real(ccs_real) :: interpol_factor
    real(ccs_real) :: dx_orth
    real(ccs_real), dimension(ndim) :: grad_phi_p
    real(ccs_real), dimension(ndim) :: grad_phi_nb
    real(ccs_real), dimension(ndim) :: x_nb, x_p, x_f, rnb_k_prime, rp_prime
    real(ccs_real), dimension(ndim) :: n

    logical :: is_boundary

    real(ccs_real) :: face_area
    real(ccs_real), dimension(ndim) :: face_norm

    real(ccs_real) :: V

    gradients(:) = 0.0_ccs_real

    call get_local_index(loc_p, index_p)
    call get_centre(loc_p, x_p)

    call count_neighbours(loc_p, nnb)
    do j = 1, nnb
      call create_face_locator(index_p, j, loc_f)
      call get_boundary_status(loc_f, is_boundary)
      call get_face_area(loc_f, face_area)
      call get_face_normal(loc_f, face_norm)

      call create_neighbour_locator(loc_p, j, loc_nb)
      call get_local_index(loc_nb, index_nb)
      if (.not. is_boundary) then
        interpol_factor = 0.5_ccs_real
        phif = interpol_factor * phi%values_ro(index_p) + (1.0_ccs_real - interpol_factor) * phi%values_ro(index_nb)

        if (phi%enable_cell_corrections) then
          call get_face_normal(loc_f, n)
          call get_centre(loc_nb, x_nb)
          call get_centre(loc_f, x_f)

          grad_phi_p = [phi%x_gradients_ro(index_p), phi%y_gradients_ro(index_p), phi%z_gradients_ro(index_p)]
          grad_phi_nb = [phi%x_gradients_ro(index_nb), phi%y_gradients_ro(index_nb), phi%z_gradients_ro(index_nb)]

          dx_orth = min(dot_product(x_f - x_p, n), dot_product(x_nb - x_f, n))
          rnb_k_prime = x_f + dx_orth * n
          rp_prime = x_f - dx_orth * n

          phif = phif + 0.5_ccs_real * (dot_product(grad_phi_nb, rnb_k_prime - x_nb) + dot_product(grad_phi_p, rp_prime - x_p))
        end if

      else
        call compute_boundary_values(phi, 0, loc_p, loc_f, face_norm, phif)
      end if

      gradients(:) = gradients(:) + phif * (face_area * face_norm(:))
    end do

    call get_volume(loc_p, V)
    gradients(:) = gradients(:) / V

  end subroutine compute_gradient_at_point

  !v Zeros the linear and fixed sources. Can be used in place of a specific implementation when
  !  there are no sources, serves as a template for any case-specific source implementation.
  !
  !  Note source routines should return "integrated" sources, i.e. SV and RV.
  module subroutine zero_sources(flow, phi, R, S)

    type(fluid), intent(in) :: flow !< Provides access to full flow field
    class(field), intent(in) :: phi !< Field being transported
    class(ccs_vector), intent(inout) :: R !< Work vector (for evaluating linear/implicit sources)
    class(ccs_vector), intent(inout) :: S !< Work vector (for evaluating fixed/explicit sources)

    real(ccs_real), dimension(:), pointer :: R_data, S_data
    integer(ccs_int) :: local_num_cells

    integer :: index_p
    type(cell_locator) :: loc_p
    real(ccs_real) :: V_p

    ! Silence warnings
    associate (foo => flow, bar => phi)
    end associate

    call get_vector_data(R, R_data)
    call get_vector_data(S, S_data)

    call get_local_num_cells(local_num_cells)
    do index_p = 1, local_num_cells
      ! XXX: Dummy implementation, use flow/phi to compute field-specific sources
      call create_cell_locator(index_p, loc_p)
      call get_volume(loc_p, V_p)
      R_data(index_p) = 0 * V_p
      S_data(index_p) = 0 * V_p
    end do

    call restore_vector_data(R, R_data)
    call restore_vector_data(S, S_data)
    call update(R)
    call update(S)

  end subroutine zero_sources

  !v Adds a fixed source term to the righthand side of the equation, expects the integrated source
  !  SV
  module subroutine add_fixed_source(S, rhs)
    use solver, only: axpy

    class(ccs_vector), intent(inout) :: S   !< The source field
    class(ccs_vector), intent(inout) :: rhs !< The righthand side vector

    call axpy(1.0_ccs_real, S, rhs)
    call update(rhs)

  end subroutine add_fixed_source

  !> Adds a linear source term to the system matrix, expects the integrated source RV
  module subroutine add_linear_source(S, M)

    use mat, only: add_matrix_diagonal

    class(ccs_vector), intent(inout) :: S !< The source field
    class(ccs_matrix), intent(inout) :: M !< The system

    call add_matrix_diagonal(S, M)

  end subroutine add_linear_source

end submodule fv_common
