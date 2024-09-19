!v Submodule file fv_common.smod
!
!  An implementation of the finite volume method
submodule(fv) fv_common
#include "ccs_macros.inc"
  use constants, only: add_mode, insert_mode
  use types, only: vector_values, matrix_values_spec, matrix_values, neighbour_locator, bc_profile, field
  use vec, only: get_vector_data, restore_vector_data, &
                 get_vector_data_readonly, restore_vector_data_readonly, &
                 create_vector_values

  use mat, only: create_matrix_values, set_matrix_values_spec_nrows, set_matrix_values_spec_ncols
  use utils, only: clear_entries, set_entry, set_row, set_col, set_values, set_mode, update
  use utils, only: debug_print, exit_print, str
  use meshing, only: count_neighbours, get_boundary_status, create_neighbour_locator, &
                     get_local_index, get_global_index, get_volume, get_distance, &
                     create_face_locator, get_face_area, get_face_normal, create_cell_locator, &
                     get_local_num_cells, get_face_interpolation, &
                     get_max_faces, get_centre
  use boundary_conditions, only: get_bc_index
  use timers, only: timer_register_start, timer_stop
  use bc_constants
  use error_codes

  implicit none

contains

  !> Computes fluxes and assign to matrix and RHS
  module subroutine compute_fluxes(phi, mf, viscosity, density, component, M, vec)
    class(field), intent(inout) :: phi
    class(field), intent(inout) :: mf
    class(field), intent(inout) :: viscosity
    class(field), intent(inout) :: density
    integer(ccs_int), intent(in) :: component
    class(ccs_matrix), intent(inout) :: M
    class(ccs_vector), intent(inout) :: vec

    integer(ccs_int) :: max_faces
    integer(ccs_int) :: n_int_cells
    real(ccs_real), dimension(:), pointer :: mf_data, viscosity_data, density_data

    associate (mf_values => mf%values)
      call dprint("CF: get mf")
      call get_vector_data(mf_values, mf_data)
      call get_vector_data(viscosity%values, viscosity_data)
      call get_vector_data(density%values, density_data)

      ! Loop over cells computing advection and diffusion fluxes
      call get_max_faces(max_faces)
      n_int_cells = max_faces + 1 ! 1 neighbour per face + central cell
      call dprint("CF: compute coeffs")
      call compute_coeffs(phi, mf_data, viscosity_data, density_data, component, n_int_cells, M, vec)

      call dprint("CF: restore mf")
      call restore_vector_data(mf_values, mf_data)
      call restore_vector_data(viscosity%values, viscosity_data)
      call restore_vector_data(density%values, density_data)
    end associate


  end subroutine compute_fluxes

  !> Computes the matrix coefficient for cells in the interior of the mesh
  subroutine compute_coeffs(phi, mf, visc, dens, component, n_int_cells, M, b)
    class(field), intent(inout) :: phi                !< scalar field structure
    real(ccs_real), dimension(:), intent(in) :: mf !< mass flux array defined at faces
    real(ccs_real), dimension(:), intent(in) :: visc !< viscosity
    real(ccs_real), dimension(:), intent(in) :: dens !< density
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
    real(ccs_real) :: adv_coeff_total, adv_coeffaF, adv_coeffaP
    real(ccs_real), dimension(ndim) :: face_normal
    real(ccs_real), dimension(ndim) :: grad_phi_p 
    real(ccs_real), dimension(ndim) :: grad_phi_nb
    real(ccs_real), dimension(ndim) :: x_nb, x_p, x_f, x_nb_prime, x_p_prime
    real(ccs_real), dimension(ndim) :: n
    real(ccs_real) :: face_value, face_correction_only

    logical :: is_boundary

    integer(ccs_int) :: index_f

    real(ccs_real) :: sgn ! Sign indicating face orientation
    real(ccs_real) :: aP, aF, bP, aPb
    real(ccs_real) :: hoe ! High-order explicit flux
    real(ccs_real) :: loe ! Low-order explicit flux
    real(ccs_real) :: dx_orth ! distance between cell centers projected to the face othogonal (used for corrections)
    real(ccs_real) :: SchmidtNo

    call set_matrix_values_spec_nrows(1_ccs_int, mat_val_spec)
    call set_matrix_values_spec_ncols(n_int_cells, mat_val_spec)
    call create_matrix_values(mat_val_spec, mat_coeffs)
    call set_mode(add_mode, mat_coeffs)

    call create_vector_values(n_int_cells, b_coeffs)
    call set_mode(add_mode, b_coeffs)

    call get_local_num_cells(local_num_cells) 
    do index_p = 1, local_num_cells
      call clear_entries(mat_coeffs)
      call clear_entries(b_coeffs)

      ! Calculate contribution from neighbours
      call create_cell_locator(index_p, loc_p)
      call get_global_index(loc_p, global_index_p)
      call count_neighbours(loc_p, nnb)

      call set_row(global_index_p, mat_coeffs)
      call set_row(global_index_p, b_coeffs)

      adv_coeff_total = 0.0_ccs_real
      diff_coeff_total = 0.0_ccs_real
      
      do j = 1, nnb
        call create_neighbour_locator(loc_p, j, loc_nb)
        call get_boundary_status(loc_nb, is_boundary)
        call create_face_locator(index_p, j, loc_f)
        call get_face_normal(loc_f, face_normal)

        call get_local_index(loc_nb, index_nb)

        call get_face_area(loc_f, face_area)
        call get_local_index(loc_f, index_f)
        SchmidtNo = phi%Schmidt
        

        if (.not. is_boundary) then
          call calc_diffusion_coeff(index_p, j, phi%enable_cell_corrections, visc(index_p), visc(index_nb), dens(index_p), dens(index_nb), SchmidtNo, diff_coeff)

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
            call calc_advection_coeff(phi, loc_f, sgn * mf(index_f), 0, adv_coeffaP, adv_coeffaF)
          type is (upwind_field)
            call calc_advection_coeff(phi, loc_f, sgn * mf(index_f), 0, adv_coeffaP, adv_coeffaF)
          type is (gamma_field)
            call calc_advection_coeff(phi, loc_f, sgn * mf(index_f), 0, loc_p, loc_nb, adv_coeffaP, adv_coeffaF)
          type is (linear_upwind_field)
            call calc_advection_coeff(phi, loc_f, sgn * mf(index_f), 0, loc_p, loc_nb, adv_coeffaP, adv_coeffaF)
          class default
            call error_abort("Invalid velocity field discretisation.")
          end select


          ! XXX: we are relying on div(u)=0 => a_P = -sum_nb a_nb
          ! adv_coeff = adv_coeff * (sgn * mf(index_f) * face_area)

          ! High-order, explicit
          aF = adv_coeffaF
          aP = adv_coeffaP
          aP = (sgn * mf(index_f) * face_area) * aP
          aF = (sgn * mf(index_f) * face_area) * aF
          aP = aP - sgn * mf(index_f) * face_area
          hoe = aP * phi%values_ro(index_p) + aF * phi%values_ro(index_nb)

          ! Excentricity correction (convective term) (Ferziger & Peric 4th ed, sec 9.7.1)
          call interpolate_field_to_face(phi, loc_f, face_value, face_correction_only)
          hoe = hoe + face_correction_only * (sgn * mf(index_f) * face_area)

          if (phi%enable_cell_corrections) then
            call get_face_normal(loc_f, n)
            call get_centre(loc_p, x_p)
            call get_centre(loc_nb, x_nb)
            call get_centre(loc_f, x_f)

            dx_orth = min(dot_product(x_f - x_p, n), dot_product(x_nb - x_f, n))
            x_nb_prime = x_f + dx_orth*n
            x_p_prime = x_f - dx_orth*n


            grad_phi_p = [ phi%x_gradients_ro(index_p), phi%y_gradients_ro(index_p), phi%z_gradients_ro(index_p) ]
            grad_phi_nb = [ phi%x_gradients_ro(index_nb), phi%y_gradients_ro(index_nb), phi%z_gradients_ro(index_nb) ]

            ! call get_face_interpolation(loc_f, interpol_factor)
            ! x_f_prime = interpol_factor * x_p + (1.0_ccs_real - interpol_factor) * x_nb
            ! grad_phi_k_prime = interpol_factor * grad_phi_p + (1.0_ccs_real - interpol_factor) * grad_phi_nb
            !hoe = hoe + dot_product(grad_phi_k_prime, x_f - x_f_prime) * (sgn * mf(index_f) * face_area)
            !hoe = hoe + 0.5_ccs_real * (dot_product(grad_phi_p, x_p_prime - x_p) + dot_product(grad_phi_nb, x_nb_prime - x_nb)) * (sgn * mf(index_f) * face_area)


            ! Non-orthogonality correction (diffusive flux) (Ferziger & Peric 4th ed, sec 9.7.2)
            hoe = hoe + diff_coeff * (dot_product(grad_phi_nb, x_nb_prime - x_nb) - dot_product(grad_phi_p, x_p_prime - x_p))
          end if

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

          loe = aP * phi%values_ro(index_p) + aF * phi%values_ro(index_nb)
          call set_entry(loe, b_coeffs) ! Explicit low-order term

          call get_global_index(loc_nb, global_index_nb)
          call set_col(global_index_nb, mat_coeffs)
          call set_entry(aF + diff_coeff, mat_coeffs)

          adv_coeff_total = adv_coeff_total + aP
          diff_coeff_total = diff_coeff_total - diff_coeff
        else
          call compute_boundary_coeffs(phi, component, loc_p, loc_f, face_normal, aPb, bP)

          call calc_diffusion_coeff(index_p, j, .false., visc(index_p), 0.0_ccs_real, dens(index_p), 0.0_ccs_real, SchmidtNo, diff_coeff)
          ! Correct boundary face distance to distance to immaginary boundary "node"
          diff_coeff = diff_coeff / 2.0_ccs_real
          
          select type (phi)
          type is (central_field)
            call calc_advection_coeff(phi, loc_f, mf(index_f), index_nb, adv_coeffaP, adv_coeffaF)
          type is (upwind_field)
            call calc_advection_coeff(phi, loc_f, mf(index_f), index_nb, adv_coeffaP, adv_coeffaF)
          type is (gamma_field)
            call calc_advection_coeff(phi, loc_f, mf(index_f), index_nb, loc_p, loc_nb, adv_coeffaP, adv_coeffaF)
          type is  (linear_upwind_field)
            call calc_advection_coeff(phi, loc_f, mf(index_f), index_nb, loc_p, loc_nb, adv_coeffaP, adv_coeffaF)
          class default
            call error_abort("Invalid velocity field discretisation.")
          end select
          aF = adv_coeffaF
          aP = adv_coeffaP
          aP = aP * (mf(index_f) * face_area)
          aF = aF * (mf(index_f) * face_area)
          aP = aP - mf(index_f) * face_area
          call set_entry(-(aP * phi%values_ro(index_p) + aF * (aPb * phi%values_ro(index_p) + bP)), b_coeffs)
          if (mf(index_f) > 0.0_ccs_real) then
            aP = mf(index_f) * face_area
            aF = 0.0_ccs_real
          else
            aP = 0.0_ccs_real
            aF = mf(index_f) * face_area
          end if
          aP = aP - mf(index_f) * face_area

          call set_entry(aP * phi%values_ro(index_p) + aF * (aPb * phi%values_ro(index_p) + bP), b_coeffs)

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
  pure module subroutine compute_boundary_values(phi, component, loc_p, loc_f, normal, bc_value)

    class(field), intent(inout) :: phi                      !< the field for which boundary values are being computed
    integer(ccs_int), intent(in) :: component               !< integer indicating direction of velocity field component
    type(cell_locator), intent(in) :: loc_p                 !< location of cell
    type(face_locator), intent(in) :: loc_f                 !< location of face
    real(ccs_real), dimension(ndim), intent(in) :: normal   !< boundary face normal direction
    real(ccs_real), intent(out) :: bc_value                 !< the boundary value

    real(ccs_real) :: a !< The diagonal coeff (implicit component)
    real(ccs_real) :: b !< The RHS value (explicit component)
    integer(ccs_int) :: index_p

    call compute_boundary_coeffs(phi, component, loc_p, loc_f, normal, a, b)

    call get_local_index(loc_p, index_p)
    bc_value = 0.5_ccs_real * (phi%values_ro(index_p) + (b + a * phi%values_ro(index_p)))

  end subroutine compute_boundary_values

  !> Compute the coefficients of the boundary condition
  pure module subroutine compute_boundary_coeffs(phi, component, loc_p, loc_f, normal, a, b)

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
    case (bc_type_sym)  ! XXX: Make sure this works as intended for symmetric BC.
      select case (component)
      case (0)
        parallel_component_map = [1, 1, 1]
      case (1)
        parallel_component_map = [0, 1, 1]
      case (2)
        parallel_component_map = [1, 0, 1]
      case (3)
        parallel_component_map = [1, 1, 0]
      case default
        error stop invalid_component ! Invalid component provided
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

    do i=1, n-1
      if (r .lt. profile%coordinates(i+1)) then
        coeff = (r - profile%coordinates(i)) / (profile%coordinates(i+1) - profile%coordinates(i))
        bc_value = (1-coeff) * profile%values(i) + coeff * profile%values(i+1)
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
    real(ccs_real), parameter :: density = 1.0_ccs_real 
    real(ccs_real) :: interpolation_factor

    call create_face_locator(index_p, index_nb, loc_f)
    call get_face_area(loc_f, face_area)
    call get_boundary_status(loc_f, is_boundary)

    call create_cell_locator(index_p, loc_p)
    call get_face_interpolation(loc_f, interpolation_factor)
    if (.not. is_boundary) then
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
    
    coeff = -face_area * diffusion_factor / dxmag
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
        x_nb_prime = x_f + dx_orth*n
        x_p_prime = x_f - dx_orth*n

        grad_phi_p = [ phi%x_gradients_ro(index_p), phi%y_gradients_ro(index_p), phi%z_gradients_ro(index_p) ]
        grad_phi_nb = [ phi%x_gradients_ro(index_nb), phi%y_gradients_ro(index_nb), phi%z_gradients_ro(index_nb) ]

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

        grad_phi_p = [ dpdx(index_p), dpdy(index_p), dpdz(index_p) ]
        grad_phi_nb = [ dpdx(index_nb), dpdy(index_nb), dpdz(index_nb) ]

        if (enable_cell_corrections) then
          ! Cell excentricity/non-orthogonality corrections (Ferziger & Peric 4th ed, sec 9.8, p317, eq9.67 and 9.66)
          call get_face_normal(loc_f, n)
          call get_centre(loc_p, x_p)
          call get_centre(loc_nb, x_nb)
          call get_centre(loc_f, x_f)

          dx_orth = min(dot_product(x_f - x_p, n), dot_product(x_nb - x_f, n))
          rnb_k_prime = x_f + dx_orth*n
          rp_prime = x_f - dx_orth*n
          !dx = rnb_k_prime - rp_prime 
          !dxmag = norm2(dx)
          dxmag = 2.0_ccs_real * dx_orth

          ! eq 9.66
          flux_corr = (p(index_nb) - p(index_p)) + (dot_product(grad_phi_nb, rnb_k_prime - x_nb) - dot_product(grad_phi_p, rp_prime - x_p))

          ! interpolated pressure gradient at the face (i.e.)
          flux_corr = flux_corr - 0.5_ccs_real * dot_product((grad_phi_p + grad_phi_nb), x_nb - x_p)

          flux_corr = - flux_corr / dxmag
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
  module subroutine update_gradient(phi)

    use meshing, only: get_local_num_cells

    class(field), intent(inout) :: phi !< the field whose gradients we want to update
    real(ccs_real), dimension(:), allocatable :: x_gradients
    real(ccs_real), dimension(:), allocatable :: y_gradients
    real(ccs_real), dimension(:), allocatable :: z_gradients

    real(ccs_real), dimension(:), pointer :: gradients_data    ! Data array for gradients
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: timer_index

    call timer_register_start("Compute gradient", timer_index)

    call get_local_num_cells(local_num_cells)
    allocate(x_gradients(local_num_cells))
    allocate(y_gradients(local_num_cells))
    allocate(z_gradients(local_num_cells))

    call dprint("Compute x gradient")
    call update_gradient_component(1, phi, x_gradients)
    call dprint("Compute y gradient")
    call update_gradient_component(2, phi, y_gradients)
    call dprint("Compute z gradient")
    call update_gradient_component(3, phi, z_gradients)

    call get_vector_data(phi%x_gradients, gradients_data)
    gradients_data(1:local_num_cells) = x_gradients(:)
    call restore_vector_data(phi%x_gradients, gradients_data)
    call update(phi%x_gradients) ! XXX: opportunity to overlap update with later compute (begin/compute/end)

    call get_vector_data(phi%y_gradients, gradients_data)
    gradients_data(1:local_num_cells) = y_gradients(:)
    call restore_vector_data(phi%y_gradients, gradients_data)
    call update(phi%y_gradients) ! yyy: opportunity to overlap update with later compute (begin/compute/end)

    call get_vector_data(phi%z_gradients, gradients_data)
    gradients_data(1:local_num_cells) = z_gradients(:)
    call restore_vector_data(phi%z_gradients, gradients_data)
    call update(phi%z_gradients) ! zzz: opportunity to overlap update with later compute (begin/compute/end)

    call timer_stop(timer_index)

  end subroutine update_gradient

  !> Helper subroutine to calculate a gradient component at a time.
  pure subroutine update_gradient_component(component, phi, gradients)

    integer(ccs_int), intent(in) :: component !< which vector component (i.e. direction) to update?
    class(field), intent(inout) :: phi !< the field whose gradient we want to compute
    real(ccs_real), dimension(:), intent(inout) :: gradients !< a cell-centred array of the gradient

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
    real(ccs_real) :: dx_orth
    real(ccs_real), dimension(ndim) :: grad_phi_p 
    real(ccs_real), dimension(ndim) :: grad_phi_nb
    real(ccs_real), dimension(ndim) :: x_nb, x_p, x_f, rnb_k_prime, rp_prime
    real(ccs_real), dimension(ndim) :: n

    logical :: is_boundary

    real(ccs_real) :: face_area
    real(ccs_real), dimension(ndim) :: face_norm

    real(ccs_real) :: V

    call get_local_num_cells(local_num_cells)
    do index_p = 1, local_num_cells
      grad = 0.0_ccs_int

      call create_cell_locator(index_p, loc_p)
      call count_neighbours(loc_p, nnb)
      call get_centre(loc_p, x_p)
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

            grad_phi_p = [ phi%x_gradients_ro(index_p), phi%y_gradients_ro(index_p), phi%z_gradients_ro(index_p) ]
            grad_phi_nb = [ phi%x_gradients_ro(index_nb), phi%y_gradients_ro(index_nb), phi%z_gradients_ro(index_nb) ]

            dx_orth = min(dot_product(x_f - x_p, n), dot_product(x_nb - x_f, n))
            rnb_k_prime = x_f + dx_orth*n
            rp_prime = x_f - dx_orth*n

            phif = phif + 0.5_ccs_real*(dot_product(grad_phi_nb, rnb_k_prime - x_nb) + dot_product(grad_phi_p, rp_prime - x_p))
          end if

        else
          call compute_boundary_values(phi, component, loc_p, loc_f, face_norm, phif)
        end if

        grad = grad + phif * (face_area * face_norm(component))
      end do

      call get_volume(loc_p, V)
      grad = grad / V
      gradients(index_p) = grad
    end do

  end subroutine update_gradient_component

  !v Zeros the linear and fixed sources. Can be used in place of a specific implementation when
  !  there are no sources, serves as a template for any case-specific source implementation.
  !
  !  Note source routines should return "integrated" sources, i.e. SV and RV.
  module subroutine zero_sources(flow, phi, R, S)

    type(fluid), intent(in) :: flow !< Provides access to full flow field
    class(field), intent(in) :: phi !< Field being transported
    class(ccs_vector), intent(inout) :: R !< Work vector (for evaluating linear/implicit sources)
    class(ccs_vector), intent(inout) :: S !< Work vector (for evaluating fixed/explicit sources)

    real(ccs_real), dimension(:), pointer:: R_data, S_data
    integer(ccs_int) :: local_num_cells

    integer :: index_p
    type(cell_locator) :: loc_p
    real(ccs_real) :: V_p

    ! Silence warnings
    associate(foo => flow, bar => phi)
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
    
    ! Ensure RHS is in "clean" state
    call update(rhs)
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
