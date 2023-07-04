!v Submodule file fv_common.smod
!
!  An implementation of the finite volume method

submodule(fv) fv_common
#include "ccs_macros.inc"
  use constants, only: add_mode, insert_mode
  use types, only: vector_values, matrix_values_spec, matrix_values, neighbour_locator, bc_profile
  use vec, only: get_vector_data, restore_vector_data, create_vector_values

  use mat, only: create_matrix_values, set_matrix_values_spec_nrows, set_matrix_values_spec_ncols
  use utils, only: clear_entries, set_entry, set_row, set_col, set_values, set_mode, update
  use utils, only: debug_print, exit_print, str
  use meshing, only: count_neighbours, get_boundary_status, create_neighbour_locator, &
                     get_local_index, get_global_index, get_volume, get_distance, &
                     create_face_locator, get_face_area, get_face_normal, create_cell_locator, &
                     get_local_num_cells, get_face_interpolation, &
                     get_max_faces, get_centre
  use boundary_conditions, only: get_bc_index
  use bc_constants

  implicit none

contains

  !> Computes fluxes and assign to matrix and RHS
  module subroutine compute_fluxes(phi, mf, mesh, component, M, vec)
    class(field), intent(inout) :: phi
    class(field), intent(inout) :: mf
    type(ccs_mesh), intent(in) :: mesh
    integer(ccs_int), intent(in) :: component
    class(ccs_matrix), intent(inout) :: M
    class(ccs_vector), intent(inout) :: vec

    integer(ccs_int) :: max_faces
    integer(ccs_int) :: n_int_cells
    real(ccs_real), dimension(:), pointer :: mf_data

    associate (mf_values => mf%values)
      call dprint("CF: get mf")
      call get_vector_data(mf_values, mf_data)

      ! Loop over cells computing advection and diffusion fluxes
      call get_max_faces(mesh, max_faces)
      n_int_cells = max_faces + 1 ! 1 neighbour per face + central cell
      call dprint("CF: compute coeffs")
      call compute_coeffs(phi, mf_data, mesh, component, n_int_cells, M, vec)

      call dprint("CF: restore mf")
      call restore_vector_data(mf_values, mf_data)
    end associate

  end subroutine compute_fluxes

  !> Computes the matrix coefficient for cells in the interior of the mesh
  subroutine compute_coeffs(phi, mf, mesh, component, n_int_cells, M, b)
    class(field), intent(inout) :: phi                !< scalar field structure
    real(ccs_real), dimension(:), intent(in) :: mf !< mass flux array defined at faces
    type(ccs_mesh), intent(in) :: mesh             !< Mesh structure
    integer(ccs_int), intent(in) :: component      !< integer indicating direction of velocity field component
    integer(ccs_int), intent(in) :: n_int_cells    !< number of cells in the interior of the mesh
    class(ccs_matrix), intent(inout) :: M          !< equation system matrix
    class(ccs_vector), intent(inout) :: b          !< RHS vector

    type(matrix_values_spec) :: mat_val_spec
    type(matrix_values) :: mat_coeffs
    type(vector_values) :: b_coeffs
    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb
    type(face_locator) :: loc_f
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: global_index_p, global_index_nb, index_p, index_nb
    integer(ccs_int) :: j
    integer(ccs_int) :: nnb
    real(ccs_real) :: face_area
    real(ccs_real) :: diff_coeff, diff_coeff_total
    real(ccs_real) :: adv_coeff, adv_coeff_total
    real(ccs_real), dimension(ndim) :: face_normal
    logical :: is_boundary

    integer(ccs_int) :: index_f

    real(ccs_real) :: sgn ! Sign indicating face orientation
    real(ccs_real) :: aP, aF, bP, aPb
    real(ccs_real), dimension(:), pointer :: phi_data
    real(ccs_real) :: hoe ! High-order explicit flux
    real(ccs_real) :: loe ! Low-order explicit flux

    call set_matrix_values_spec_nrows(1_ccs_int, mat_val_spec)
    call set_matrix_values_spec_ncols(n_int_cells, mat_val_spec)
    call create_matrix_values(mat_val_spec, mat_coeffs)
    call set_mode(add_mode, mat_coeffs)

    call create_vector_values(n_int_cells, b_coeffs)
    call set_mode(add_mode, b_coeffs)

    call get_local_num_cells(mesh, local_num_cells)
    do index_p = 1, local_num_cells
      call clear_entries(mat_coeffs)
      call clear_entries(b_coeffs)

      ! Calculate contribution from neighbours
      call create_cell_locator(mesh, index_p, loc_p)
      call get_global_index(loc_p, global_index_p)
      call count_neighbours(loc_p, nnb)

      call set_row(global_index_p, mat_coeffs)
      call set_row(global_index_p, b_coeffs)

      adv_coeff_total = 0.0_ccs_real
      diff_coeff_total = 0.0_ccs_real

      do j = 1, nnb
        call create_neighbour_locator(loc_p, j, loc_nb)
        call get_boundary_status(loc_nb, is_boundary)
        call create_face_locator(mesh, index_p, j, loc_f)
        call get_face_normal(loc_f, face_normal)

        call get_local_index(loc_nb, index_nb)

        call get_face_area(loc_f, face_area)
        call get_local_index(loc_f, index_f)

        if (.not. is_boundary) then
          diff_coeff = calc_diffusion_coeff(index_p, j, mesh)

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
            call calc_advection_coeff(phi, loc_f, sgn * mf(index_f), 0, adv_coeff)
          type is (upwind_field)
            call calc_advection_coeff(phi, loc_f, sgn * mf(index_f), 0, adv_coeff)
          class default
            call error_abort("Invalid velocity field discretisation.")
          end select

          call get_vector_data(phi%values, phi_data)

          ! XXX: we are relying on div(u)=0 => a_P = -sum_nb a_nb
          ! adv_coeff = adv_coeff * (sgn * mf(index_f) * face_area)

          ! High-order, explicit
          aF = adv_coeff
          aP = 1.0_ccs_real - aF
          aP = (sgn * mf(index_f) * face_area) * aP
          aF = (sgn * mf(index_f) * face_area) * aF
          aP = aP - sgn * mf(index_f) * face_area
          hoe = aP * phi_data(index_p) + aF * phi_data(index_nb)
          call set_entry(-hoe, b_coeffs)

          ! Low-order
          if ((sgn * mf(index_f)) > 0.0_ccs_real) then
            aP = sgn * mf(index_f) * face_area
            aF = 0.0_ccs_real
          else
            aP = 0.0_ccs_real
            aF = sgn * mf(index_f) * face_area
          end if
          aP = aP - sgn * mf(index_f) * face_area

          loe = aP * phi_data(index_p) + aF * phi_data(index_nb)
          call set_entry(loe, b_coeffs) ! Explicit low-order term

          call get_global_index(loc_nb, global_index_nb)
          call set_col(global_index_nb, mat_coeffs)
          call set_entry(aF + diff_coeff, mat_coeffs)

          adv_coeff_total = adv_coeff_total + aP
          diff_coeff_total = diff_coeff_total - diff_coeff

          call restore_vector_data(phi%values, phi_data)
        else
          call compute_boundary_coeffs(phi, component, loc_p, loc_f, face_normal, aPb, bP)

          diff_coeff = calc_diffusion_coeff(index_p, j, mesh)
          ! Correct boundary face distance to distance to immaginary boundary "node"
          diff_coeff = diff_coeff / 2.0_ccs_real

          select type (phi)
          type is (central_field)
            call calc_advection_coeff(phi, loc_f, mf(index_f), index_nb, adv_coeff)
          type is (upwind_field)
            call calc_advection_coeff(phi, loc_f, mf(index_f), index_nb, adv_coeff)
          class default
            call error_abort("Invalid velocity field discretisation.")
          end select
          aF = adv_coeff
          aP = 1.0_ccs_real - aF
          aP = aP * (mf(index_f) * face_area)
          aF = aF * (mf(index_f) * face_area)
          aP = aP - mf(index_f) * face_area
          call get_vector_data(phi%values, phi_data)
          call set_entry(-(aP * phi_data(index_p) + aF * (aPb * phi_data(index_p) + bP)), b_coeffs)
          if (mf(index_f) > 0.0_ccs_real) then
            aP = mf(index_f) * face_area
            aF = 0.0_ccs_real
          else
            aP = 0.0_ccs_real
            aF = mf(index_f) * face_area
          end if
          aP = aP - mf(index_f) * face_area

          call set_entry(aP * phi_data(index_p) + aF * (aPb * phi_data(index_p) + bP), b_coeffs)
          call restore_vector_data(phi%values, phi_data)

          call set_entry(-(aF + diff_coeff) * bP, b_coeffs)

          adv_coeff_total = adv_coeff_total + aP + aPb * aF
          diff_coeff_total = diff_coeff_total - diff_coeff + aPb * diff_coeff
        end if
      end do

      call set_values(b_coeffs, b)
      call set_col(global_index_p, mat_coeffs)
      call set_entry((adv_coeff_total + diff_coeff_total), mat_coeffs)
      call set_values(mat_coeffs, M)
    end do

    deallocate (mat_coeffs%global_row_indices)
    deallocate (mat_coeffs%global_col_indices)
    deallocate (mat_coeffs%values)
  end subroutine compute_coeffs

  !> Computes the value of the scalar field on the boundary
  module subroutine compute_boundary_values(phi, component, loc_p, loc_f, normal, bc_value)

    class(field), intent(inout) :: phi                      !< the field for which boundary values are being computed
    integer(ccs_int), intent(in) :: component               !< integer indicating direction of velocity field component
    type(cell_locator), intent(in) :: loc_p                 !< location of cell
    type(face_locator), intent(in) :: loc_f                 !< location of face
    real(ccs_real), dimension(ndim), intent(in) :: normal   !< boundary face normal direction
    real(ccs_real), intent(out) :: bc_value                 !< the boundary value

    real(ccs_real) :: a !< The diagonal coeff (implicit component)
    real(ccs_real) :: b !< The RHS value (explicit component)
    integer(ccs_int) :: index_p
    real(ccs_real), dimension(:), pointer :: phi_data

    call compute_boundary_coeffs(phi, component, loc_p, loc_f, normal, a, b)

    call get_local_index(loc_p, index_p)
    call get_vector_data(phi%values, phi_data)
    bc_value = 0.5_ccs_real * (phi_data(index_p) + (b + a * phi_data(index_p)))
    call restore_vector_data(phi%values, phi_data)

  end subroutine compute_boundary_values

  !> Compute the coefficients of the boundary condition
  module subroutine compute_boundary_coeffs(phi, component, loc_p, loc_f, normal, a, b)

    class(field), intent(inout) :: phi                      !< the field for which boundary values are being computed
    integer(ccs_int), intent(in) :: component               !< integer indicating direction of velocity field component
    type(cell_locator), intent(in) :: loc_p                 !< location of cell
    type(face_locator), intent(in) :: loc_f                 !< location of face
    real(ccs_real), dimension(ndim), intent(in) :: normal   !< boundary face normal direction
    real(ccs_real), intent(out) :: a                        !< The diagonal coeff (implicit)
    real(ccs_real), intent(out) :: b                        !< The RHS entry (explicit)

    ! local variables
    integer(ccs_int) :: index_bc
    integer(ccs_int) :: index_nb
    integer(ccs_int) :: index_p
    type(neighbour_locator) :: loc_nb
    integer(ccs_int) :: i
    real(ccs_real), dimension(ndim) :: dx
    real(ccs_real), dimension(ndim) :: x
    real(ccs_real), dimension(ndim) :: parallel_component_map
    real(ccs_real), dimension(ndim) :: phi_face_parallel_component
    real(ccs_real) :: phi_face_parallel_component_norm
    real(ccs_real) :: phi_face_parallel_component_portion
    real(ccs_real) :: normal_norm
    real(ccs_real) :: dxmag
    real(ccs_real) :: bc_value
    real(ccs_real), dimension(:), pointer :: x_gradients_data    ! Data array for x gradient
    real(ccs_real), dimension(:), pointer :: y_gradients_data    ! Data array for y gradient
    real(ccs_real), dimension(:), pointer :: z_gradients_data    ! Data array for z gradient

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

      call get_vector_data(phi%x_gradients, x_gradients_data)
      call get_vector_data(phi%y_gradients, y_gradients_data)
      call get_vector_data(phi%z_gradients, z_gradients_data)

      a = 1.0_ccs_real
      b = 2.0_ccs_real * (x_gradients_data(index_p) * dx(1) + y_gradients_data(index_p) * dx(2) + z_gradients_data(index_p) * dx(3))

      call restore_vector_data(phi%x_gradients, x_gradients_data)
      call restore_vector_data(phi%y_gradients, y_gradients_data)
      call restore_vector_data(phi%z_gradients, z_gradients_data)

    case (bc_type_sym)  ! XXX: Make sure this works as intended for symmetric BC.
      select case (component)
      case (0)
        parallel_component_map = (/1, 1, 1/)
      case (1)
        parallel_component_map = (/0, 1, 1/)
      case (2)
        parallel_component_map = (/1, 0, 1/)
      case (3)
        parallel_component_map = (/1, 1, 0/)
      case default
        call error_abort("invalid component provided " // str(component))
      end select
      ! Only keep the components of phi that are parallel to the surface
      phi_face_parallel_component_norm = 0
      normal_norm = 0
      do i = 1, ndim
        phi_face_parallel_component(i) = parallel_component_map(i) * normal(i)
        phi_face_parallel_component_norm = phi_face_parallel_component_norm + &
                                           phi_face_parallel_component(i) * phi_face_parallel_component(i)
        normal_norm = normal_norm + normal(i) * normal(i)
      end do
      phi_face_parallel_component_portion = sqrt(phi_face_parallel_component_norm / normal_norm)

      ! Get value of phi at boundary cell
      a = phi_face_parallel_component_portion
      b = 0.0_ccs_real
    case (bc_type_neumann)
      call get_distance(loc_p, loc_f, dx)
      dxmag = sqrt(sum(dx**2))

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

      call error_abort("unknown bc type " // str(phi%bcs%bc_types(index_bc)))
    end select

  end subroutine

  !> Linear interpolate of BC profile 
  module subroutine get_value_from_bc_profile(x, profile, bc_value)
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

    do i=1, n-1
      if (r .lt. profile%coordinates(i+1)) then
        coeff = (r - profile%coordinates(i)) / (profile%coordinates(i+1) - profile%coordinates(i))
        bc_value = (1-coeff) * profile%values(i) + coeff * profile%values(i+1)
        return
      end if
    end do

  end subroutine

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

    call create_face_locator(mesh, index_p, index_nb, loc_f)
    call get_face_area(loc_f, face_area)
    call get_boundary_status(loc_f, is_boundary)

    call create_cell_locator(mesh, index_p, loc_p)
    if (.not. is_boundary) then
      call create_neighbour_locator(loc_p, index_nb, loc_nb)
      call get_distance(loc_p, loc_nb, dx)
    else
      call get_distance(loc_p, loc_f, dx)
    end if
    dxmag = sqrt(sum(dx**2))

    coeff = -face_area * diffusion_factor / dxmag
  end function calc_diffusion_coeff

  !> Calculates mass flux across given face. Note: assumes rho = 1 and uniform grid
  module function calc_mass_flux_uvw(u_field, v_field, w_field, p, dpdx, dpdy, dpdz, invAu, invAv, invAw, loc_f) result(flux)
    class(field), intent(inout) :: u_field
    class(field), intent(inout) :: v_field
    class(field), intent(inout) :: w_field
    real(ccs_real), dimension(:), intent(in) :: p                   !< array containing pressure
    real(ccs_real), dimension(:), intent(in) :: dpdx, dpdy, dpdz    !< arrays containing pressure gradient in x, y and z
    real(ccs_real), dimension(:), intent(in) :: invAu, invAv, invAw !< arrays containing inverse momentum diagonal in x, y and z
    type(face_locator), intent(in) :: loc_f                         !< face locator

    real(ccs_real) :: flux                                          !< The flux across the boundary

    ! Local variables
    logical :: is_boundary                         ! Boundary indicator
    type(cell_locator) :: loc_p                    ! Primary cell locator
    type(neighbour_locator) :: loc_nb              ! Neighbour cell locator
    integer(ccs_int) :: index_nb                   ! Neighbour cell index
    real(ccs_real) :: flux_corr                    ! Flux correction
    real(ccs_real), dimension(ndim) :: face_normal ! (local) face-normal array
    real(ccs_real) :: u_bc, v_bc, w_bc             ! values of u, v and w at boundary
    real(ccs_real), dimension(:), pointer :: u_data, v_data, w_data ! the vector data for u, v and w
    integer(ccs_int), parameter :: x_direction = 1, y_direction = 2, z_direction = 3
    real(ccs_real) :: interpol_factor

    call get_boundary_status(loc_f, is_boundary)

    associate (mesh => loc_f%mesh, &
               index_p => loc_f%index_p, &
               j => loc_f%cell_face_ctr)

      call create_cell_locator(mesh, index_p, loc_p)
      call create_neighbour_locator(loc_p, j, loc_nb)
      call get_local_index(loc_nb, index_nb)

      call get_face_normal(loc_f, face_normal)
      if (.not. is_boundary) then

        ! XXX: this is likely expensive inside a loop...
        call get_vector_data(u_field%values, u_data)
        call get_vector_data(v_field%values, v_data)
        call get_vector_data(w_field%values, w_data)

        call get_face_interpolation(loc_f, interpol_factor)

        flux = ((interpol_factor * u_data(index_p) + (1.0_ccs_real - interpol_factor) * u_data(index_nb)) * face_normal(x_direction) &
              + (interpol_factor * v_data(index_p) + (1.0_ccs_real - interpol_factor) * v_data(index_nb)) * face_normal(y_direction) &
              + (interpol_factor * w_data(index_p) + (1.0_ccs_real - interpol_factor) * w_data(index_nb)) * face_normal(z_direction))

        call restore_vector_data(u_field%values, u_data)
        call restore_vector_data(v_field%values, v_data)
        call restore_vector_data(w_field%values, w_data)

        if (index_p > index_nb) then
          ! XXX: making convention to point from low to high cell.
          flux = -flux
        end if

        flux_corr = calc_mass_flux(p, dpdx, dpdy, dpdz, invAu, invAv, invAw, loc_f)
        flux = flux + flux_corr
      else
        call compute_boundary_values(u_field, x_direction, loc_p, loc_f, face_normal, u_bc)
        call compute_boundary_values(v_field, y_direction, loc_p, loc_f, face_normal, v_bc)
        call compute_boundary_values(w_field, z_direction, loc_p, loc_f, face_normal, w_bc)
        flux = u_bc * face_normal(x_direction) + v_bc * face_normal(y_direction) + w_bc * face_normal(z_direction)
      end if
    end associate

  end function calc_mass_flux_uvw

  ! Computes Rhie-Chow correction
  module function calc_mass_flux_no_uvw(p, dpdx, dpdy, dpdz, invAu, invAv, invAw, loc_f) result(flux)
    real(ccs_real), dimension(:), intent(in) :: p                   !< array containing pressure
    real(ccs_real), dimension(:), intent(in) :: dpdx, dpdy, dpdz    !< arrays containing pressure gradient in x, y and z
    real(ccs_real), dimension(:), intent(in) :: invAu, invAv, invAw !< arrays containing inverse momentum diagonal in x, y and z
    type(face_locator), intent(in) :: loc_f                         !< face locator

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

    real(ccs_real) :: uSwitch, vSwitch, wSwitch
    real(ccs_real) :: problem_dim
    real(ccs_real) :: interpol_factor

    call get_boundary_status(loc_f, is_boundary)

    associate (mesh => loc_f%mesh, &
               index_p => loc_f%index_p, &
               j => loc_f%cell_face_ctr)

      call create_cell_locator(mesh, index_p, loc_p)
      call create_neighbour_locator(loc_p, j, loc_nb)
      call get_local_index(loc_nb, index_nb)

      call get_face_normal(loc_f, face_normal)
      if (.not. is_boundary) then

        !
        ! Rhie-Chow correction from Ferziger & Peric
        !
        call get_face_interpolation(loc_f, interpol_factor)
        call get_distance(loc_p, loc_nb, dx)
        dxmag = sqrt(sum(dx**2))
        call get_face_normal(loc_f, face_normal)
        flux_corr = -(p(index_nb) - p(index_p)) / dxmag

        flux_corr = flux_corr + ((interpol_factor * dpdx(index_p) + (1.0_ccs_real - interpol_factor) * dpdx(index_nb)) * face_normal(x_direction) &
                + (interpol_factor * dpdy(index_p) + (1.0_ccs_real - interpol_factor) * dpdy(index_nb)) * face_normal(y_direction) &
                 + (interpol_factor * dpdz(index_p) + (1.0_ccs_real - interpol_factor) * dpdz(index_nb)) * face_normal(z_direction))

        call get_volume(loc_p, Vp)
        call get_volume(loc_nb, V_nb)
        Vf = interpol_factor * Vp + (1.0_ccs_real - interpol_factor) * V_nb

        ! This is probably not quite right ...
        uSwitch = 1.0_ccs_real
        vSwitch = 1.0_ccs_real
        wSwitch = 1.0_ccs_real
        problem_dim = uSwitch + vSwitch + wSwitch
        invAp = (uSwitch * invAu(index_p) + vSwitch * invAv(index_p) + wSwitch * invAw(index_p)) / problem_dim
        invA_nb = (uSwitch * invAu(index_nb) + vSwitch * invAv(index_nb) + wSwitch * invAw(index_nb)) / problem_dim
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

    use meshing, only: get_total_num_cells

    type(ccs_mesh), intent(in) :: mesh !< the mesh
    class(field), intent(inout) :: phi !< the field whose gradients we want to update

    call dprint("Compute x gradient")
    call update_gradient_component(mesh, 1, phi, phi%x_gradients)
    call update(phi%x_gradients) ! XXX: opportunity to overlap update with later compute (begin/compute/end)
    call dprint("Compute y gradient")
    call update_gradient_component(mesh, 2, phi, phi%y_gradients)
    call update(phi%y_gradients) ! XXX: opportunity to overlap update with later compute (begin/compute/end)
    call dprint("Compute z gradient")
    call update_gradient_component(mesh, 3, phi, phi%z_gradients)
    call update(phi%z_gradients) ! XXX: opportunity to overlap update with later compute (begin/compute/end)

  end subroutine update_gradient

  !> Helper subroutine to calculate a gradient component at a time.
  subroutine update_gradient_component(mesh, component, phi, gradients)

    type(ccs_mesh), intent(in) :: mesh !< the mesh
    integer(ccs_int), intent(in) :: component !< which vector component (i.e. direction) to update?
    class(field), intent(inout) :: phi !< the field whose gradient we want to compute
    class(ccs_vector), intent(inout) :: gradients !< a cell-centred array of the gradient

    type(vector_values) :: grad_values
    real(ccs_real), dimension(:), pointer :: phi_data
    real(ccs_real) :: grad

    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: index_p
    integer(ccs_int) :: j
    type(cell_locator) :: loc_p
    type(face_locator) :: loc_f
    type(neighbour_locator) :: loc_nb

    integer(ccs_int) :: nnb
    integer(ccs_int) :: index_nb

    real(ccs_real) :: phif
    real(ccs_real) :: interpol_factor

    logical :: is_boundary

    real(ccs_real) :: face_area
    real(ccs_real), dimension(ndim) :: face_norm

    real(ccs_real) :: V
    integer(ccs_int) :: global_index_p

    call create_vector_values(1_ccs_int, grad_values)
    call set_mode(insert_mode, grad_values)

    call get_local_num_cells(mesh, local_num_cells)
    do index_p = 1, local_num_cells
      call clear_entries(grad_values)

      grad = 0.0_ccs_int

      call create_cell_locator(mesh, index_p, loc_p)
      call count_neighbours(loc_p, nnb)
      do j = 1, nnb
        call create_face_locator(mesh, index_p, j, loc_f)
        call get_boundary_status(loc_f, is_boundary)
        call get_face_area(loc_f, face_area)
        call get_face_normal(loc_f, face_norm)

        call create_neighbour_locator(loc_p, j, loc_nb)
        call get_local_index(loc_nb, index_nb)
        if (.not. is_boundary) then
          call get_vector_data(phi%values, phi_data)
          call get_face_interpolation(loc_f, interpol_factor)
          phif = interpol_factor * phi_data(index_p) + (1.0_ccs_real - interpol_factor) * phi_data(index_nb)
          call restore_vector_data(phi%values, phi_data)
        else
          call compute_boundary_values(phi, component, loc_p, loc_f, face_norm, phif)
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

    deallocate (grad_values%global_indices)
    deallocate (grad_values%values)

  end subroutine update_gradient_component

  !> Adds a fixed source term to the righthand side of the equation
  module subroutine add_fixed_source(mesh, S, rhs)

    type(ccs_mesh), intent(in) :: mesh     !< The mesh
    class(ccs_vector), intent(inout) :: S      !< The source field
    class(ccs_vector), intent(inout) :: rhs !< The righthand side vector

    real(ccs_real), dimension(:), pointer :: S_data
    real(ccs_real), dimension(:), pointer :: rhs_data

    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: index_p
    type(cell_locator) :: loc_p
    real(ccs_real) :: V_p
    
    call get_vector_data(S, S_data)
    call get_vector_data(rhs, rhs_data)
    call get_local_num_cells(mesh, local_num_cells)
    do index_p = 1, local_num_cells
       call create_cell_locator(mesh, index_p, loc_p)
       call get_volume(loc_p, V_p)

       rhs_data(index_p) = rhs_data(index_p) + S_data(index_p) * V_p
    end do
    call restore_vector_data(S, S_data)
    call restore_vector_data(rhs, rhs_data)

    call update(rhs)
    
  end subroutine add_fixed_source

  !> Adds a linear source term to the system matrix
  module subroutine add_linear_source(mesh, S, M)

    use mat, only: add_matrix_diagonal
    
    type(ccs_mesh), intent(in) :: mesh    !< The mesh
    class(ccs_vector), intent(inout) :: S !< The source field
    class(ccs_matrix), intent(inout) :: M !< The system

    real(ccs_real), dimension(:), pointer :: S_data

    real(ccs_real), dimension(:), allocatable :: S_store
    
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: index_p
    type(cell_locator) :: loc_p
    real(ccs_real) :: V_p
    
    call get_local_num_cells(mesh, local_num_cells)
    allocate(S_store(local_num_cells))

    call get_vector_data(S, S_data)
    do index_p = 1, local_num_cells
       call create_cell_locator(mesh, index_p, loc_p)
       call get_volume(loc_p, V_p)

       S_store(index_p) = S_data(index_p)
       S_data(index_p) = S_data(index_p) * V_p
    end do
    call restore_vector_data(S, S_data)
    call update(S)

    call add_matrix_diagonal(S, M)
    
    call get_vector_data(S, S_data)
    do index_p = 1, local_num_cells
       S_data(index_p) = S_store(index_p)
    end do
    call restore_vector_data(S, S_data)
    call update(S)

    deallocate(S_store)
    
  end subroutine add_linear_source

end submodule fv_common
