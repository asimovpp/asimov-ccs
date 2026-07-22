!v Submodule file fv_equations_scalar_transport.smod
submodule(fv_equations) fv_equations_scalar_transport

  use fv, only: calc_diffusion_coeff, compute_boundary_coeffs
  use meshing, only: count_neighbours, create_face_locator, create_neighbour_locator, &
                     get_boundary_status, get_centre, get_face_area, get_face_interpolation, &
                     get_face_normal, get_global_index, get_local_index
  use types, only: face_locator, neighbour_locator

  implicit none

contains

  module subroutine scalar_transport_init(self, max_faces, component)
    class(scalar_transport_equation), intent(inout) :: self
    integer(ccs_int), intent(in) :: max_faces
    integer(ccs_int), intent(in), optional :: component

    if (present(component)) then
      self%component = component
    else
      self%component = 0_ccs_int
    end if

    call ensure_scalar_transport_capacity(self, max_faces)
    call ensure_row_capacity(self%row, max_faces + 1_ccs_int)

    self%payload%nnb = 0_ccs_int
    self%payload%index_p = 0_ccs_int
    self%payload%global_index_p = 0_ccs_int

  end subroutine scalar_transport_init

  module subroutine scalar_transport_set_advection(self, kernel)
    class(scalar_transport_equation), intent(inout) :: self
    class(advection_kernel), intent(in) :: kernel

    if (allocated(self%adv_kernel)) then
      deallocate (self%adv_kernel)
    end if
    allocate (self%adv_kernel, source = kernel)

  end subroutine scalar_transport_set_advection

  module subroutine scalar_transport_gather(self, phi, loc_p, mass_flux, viscosity, density)
    class(scalar_transport_equation), intent(inout) :: self
    class(field), intent(in) :: phi
    type(cell_locator), intent(in) :: loc_p
    real(ccs_real), dimension(:), intent(in) :: mass_flux
    real(ccs_real), dimension(:), intent(in) :: viscosity
    real(ccs_real), dimension(:), intent(in) :: density

    type(neighbour_locator) :: loc_nb
    type(face_locator) :: loc_f
    integer(ccs_int) :: index_p
    integer(ccs_int) :: global_index_p
    integer(ccs_int) :: nnb
    integer(ccs_int) :: j
    integer(ccs_int) :: index_nb
    integer(ccs_int) :: global_index_nb
    integer(ccs_int) :: index_f
    logical :: is_boundary_f
    real(ccs_real) :: phi_p
    real(ccs_real) :: phi_nb
    real(ccs_real) :: flux_f
    real(ccs_real) :: flux_mass
    real(ccs_real) :: area_f
    real(ccs_real) :: diff_coeff_f
    real(ccs_real) :: interpolation_factor_f
    real(ccs_real) :: aPb_f
    real(ccs_real) :: bP_f
    real(ccs_real) :: correction_f
    real(ccs_real) :: schmidt_no
    real(ccs_real) :: dx_orth
    real(ccs_real), dimension(ndim) :: normal_f
    real(ccs_real), dimension(ndim) :: grad_phi_p
    real(ccs_real), dimension(ndim) :: grad_phi_nb
    real(ccs_real), dimension(ndim) :: x_p
    real(ccs_real), dimension(ndim) :: x_nb
    real(ccs_real), dimension(ndim) :: x_f
    real(ccs_real), dimension(ndim) :: x_nb_prime
    real(ccs_real), dimension(ndim) :: x_p_prime

    if (.not. allocated(self%adv_kernel)) then
      error stop "scalar_transport_equation%gather: advection kernel not configured"
    end if

    call get_local_index(loc_p, index_p)
    call get_global_index(loc_p, global_index_p)
    call count_neighbours(loc_p, nnb)

    call ensure_scalar_transport_capacity(self, nnb)

    self%payload%index_p = index_p
    self%payload%global_index_p = global_index_p
    self%payload%nnb = nnb

    phi_p = phi%values_ro(index_p)
    schmidt_no = phi%Schmidt

    do j = 1, nnb
      call create_neighbour_locator(loc_p, j, loc_nb)
      call get_boundary_status(loc_nb, is_boundary_f)
      call create_face_locator(index_p, j, loc_f)
      call get_local_index(loc_f, index_f)
      call get_face_area(loc_f, area_f)

      if (.not. is_boundary_f) then
        call get_face_interpolation(loc_f, interpolation_factor_f)
        self%payload%interpolation_factor_f(j) = interpolation_factor_f
      else
        self%payload%interpolation_factor_f(j) = 0.5_ccs_real
      end if

      self%payload%is_boundary_f(j) = is_boundary_f
      self%payload%gradients_f(:, :, j) = 0.0_ccs_real
      self%payload%rvecs_f(:, :, j) = 0.0_ccs_real
      self%payload%aPb_f(j) = 0.0_ccs_real
      self%payload%bP_f(j) = 0.0_ccs_real

      flux_mass = mass_flux(index_f)
      correction_f = 0.0_ccs_real

      if (.not. is_boundary_f) then
        call get_local_index(loc_nb, index_nb)
        call get_global_index(loc_nb, global_index_nb)
        self%payload%global_indices_nb(j) = global_index_nb

        if (index_nb < index_p) flux_mass = -flux_mass

        phi_nb = phi%values_ro(index_nb)

        if (phi%enable_cell_corrections) then
          call get_face_normal(loc_f, normal_f)
          call get_centre(loc_p, x_p)
          call get_centre(loc_nb, x_nb)
          call get_centre(loc_f, x_f)

          dx_orth = min(dot_product(x_f - x_p, normal_f), dot_product(x_nb - x_f, normal_f))
          x_nb_prime = x_f + dx_orth * normal_f
          x_p_prime = x_f - dx_orth * normal_f

          grad_phi_p = [phi%x_gradients_ro(index_p), phi%y_gradients_ro(index_p), phi%z_gradients_ro(index_p)]
          grad_phi_nb = [phi%x_gradients_ro(index_nb), phi%y_gradients_ro(index_nb), phi%z_gradients_ro(index_nb)]

          self%payload%gradients_f(:, 1, j) = grad_phi_p
          self%payload%gradients_f(:, 2, j) = grad_phi_nb
          self%payload%rvecs_f(:, 1, j) = x_p_prime - x_p
          self%payload%rvecs_f(:, 2, j) = x_nb_prime - x_nb

          correction_f = 0.5_ccs_real * (dot_product(grad_phi_p, self%payload%rvecs_f(:, 1, j)) + &
                                         dot_product(grad_phi_nb, self%payload%rvecs_f(:, 2, j)))
        end if

        call calc_diffusion_coeff(index_p, j, phi%enable_cell_corrections, &
                                  viscosity(index_p), viscosity(index_nb), &
                                  density(index_p), density(index_nb), &
                                  schmidt_no, diff_coeff_f)
      else
        self%payload%global_indices_nb(j) = -1_ccs_int

        call get_face_normal(loc_f, normal_f)
        call compute_boundary_coeffs(phi, self%component, loc_p, loc_f, normal_f, aPb_f, bP_f)

        self%payload%aPb_f(j) = aPb_f
        self%payload%bP_f(j) = bP_f

        phi_nb = aPb_f * phi_p + bP_f

        call calc_diffusion_coeff(index_p, j, .false., &
                                  viscosity(index_p), 0.0_ccs_real, &
                                  density(index_p), 0.0_ccs_real, &
                                  schmidt_no, diff_coeff_f)
        diff_coeff_f = diff_coeff_f / 2.0_ccs_real
      end if

      flux_f = flux_mass * area_f
      self%payload%flux_f(j) = flux_f
      self%payload%diff_coeff_f(j) = diff_coeff_f
      self%payload%phi_values_f(1, j) = phi_p
      self%payload%phi_values_f(2, j) = phi_nb
      self%payload%correction_f(j) = correction_f
    end do

  end subroutine scalar_transport_gather

  module subroutine scalar_transport_apply(self)
    class(scalar_transport_equation), intent(inout) :: self

    real(ccs_real), dimension(2) :: diff_coeffs
    real(ccs_real), dimension(2) :: adv_coeffs
    real(ccs_real) :: diff_expl
    real(ccs_real) :: adv_expl
    real(ccs_real) :: a_P
    real(ccs_real) :: a_F
    real(ccs_real) :: rhs
    integer(ccs_int) :: j
    integer(ccs_int) :: global_row_index

    global_row_index = self%payload%global_index_p
    a_P = 0.0_ccs_real
    rhs = 0.0_ccs_real

    call self%prepare_row(global_row_index, self%payload%nnb + 1_ccs_int)

    do j = 1, self%payload%nnb
      diff_coeffs = self%diff_kernel%eval_coeffs(self%payload%diff_coeff_f(j))
      diff_expl = self%diff_kernel%eval_explicit(self%payload%phi_values_f(:, j), &
                                                 self%payload%diff_coeff_f(j), &
                                                 self%payload%interpolation_factor_f(j), &
                                                 self%payload%rvecs_f(:, :, j), &
                                                 self%payload%gradients_f(:, :, j))

      adv_coeffs = self%adv_kernel%eval_coeffs(self%payload%flux_f(j))
      adv_expl = self%adv_kernel%eval_explicit(self%payload%phi_values_f(:, j), &
                                               self%payload%flux_f(j), &
                                               self%payload%interpolation_factor_f(j), &
                                               self%payload%rvecs_f(:, :, j), &
                                               self%payload%gradients_f(:, :, j))

      rhs = rhs - diff_expl - adv_expl - self%payload%correction_f(j) * self%payload%flux_f(j)
      a_P = a_P + adv_coeffs(1) + diff_coeffs(1)
      a_F = adv_coeffs(2) + diff_coeffs(2)

      if (.not. self%payload%is_boundary_f(j)) then
        call self%add_row_entry(self%payload%global_indices_nb(j), a_F)
      else
        a_P = a_P + self%payload%aPb_f(j) * a_F
        rhs = rhs - a_F * self%payload%bP_f(j)
      end if
    end do

    call self%add_row_entry(global_row_index, a_P)
    call self%set_rhs(rhs)

  end subroutine scalar_transport_apply

  module subroutine ensure_scalar_transport_capacity(self, required)
    class(scalar_transport_equation), intent(inout) :: self
    integer(ccs_int), intent(in) :: required

    integer(ccs_int) :: allocate_size
    integer(ccs_int) :: adv_width
    integer(ccs_int) :: diff_width
    integer(ccs_int) :: stencil_width

    if (.not. allocated(self%payload)) then
      allocate (scalar_transport_payload :: self%payload)
    end if

    diff_width = self%diff_kernel%get_width()
    if (allocated(self%adv_kernel)) then
      adv_width = self%adv_kernel%get_width()
    else
      adv_width = 1_ccs_int
    end if
    stencil_width = max(diff_width, adv_width) + 1_ccs_int

    if (stencil_width /= 2_ccs_int) then
      error stop "scalar_transport_equation%ensure_capacity: kernel width > 1 not yet supported"
    end if

    allocate_size = max(1_ccs_int, required)
    if (allocate_size <= self%capacity) return

    if (allocated(self%payload%global_indices_nb)) deallocate (self%payload%global_indices_nb)
    if (allocated(self%payload%is_boundary_f)) deallocate (self%payload%is_boundary_f)
    if (allocated(self%payload%flux_f)) deallocate (self%payload%flux_f)
    if (allocated(self%payload%diff_coeff_f)) deallocate (self%payload%diff_coeff_f)
    if (allocated(self%payload%interpolation_factor_f)) deallocate (self%payload%interpolation_factor_f)
    if (allocated(self%payload%aPb_f)) deallocate (self%payload%aPb_f)
    if (allocated(self%payload%bP_f)) deallocate (self%payload%bP_f)
    if (allocated(self%payload%correction_f)) deallocate (self%payload%correction_f)
    if (allocated(self%payload%phi_values_f)) deallocate (self%payload%phi_values_f)
    if (allocated(self%payload%gradients_f)) deallocate (self%payload%gradients_f)
    if (allocated(self%payload%rvecs_f)) deallocate (self%payload%rvecs_f)

    allocate (self%payload%global_indices_nb(allocate_size))
    allocate (self%payload%is_boundary_f(allocate_size))
    allocate (self%payload%flux_f(allocate_size))
    allocate (self%payload%diff_coeff_f(allocate_size))
    allocate (self%payload%interpolation_factor_f(allocate_size))
    allocate (self%payload%aPb_f(allocate_size))
    allocate (self%payload%bP_f(allocate_size))
    allocate (self%payload%correction_f(allocate_size))
    allocate (self%payload%phi_values_f(stencil_width, allocate_size))
    allocate (self%payload%gradients_f(ndim, stencil_width, allocate_size))
    allocate (self%payload%rvecs_f(ndim, stencil_width, allocate_size))

    self%capacity = allocate_size

    self%payload%global_indices_nb = 0_ccs_int
    self%payload%is_boundary_f = .false.
    self%payload%flux_f = 0.0_ccs_real
    self%payload%diff_coeff_f = 0.0_ccs_real
    self%payload%interpolation_factor_f = 0.0_ccs_real
    self%payload%aPb_f = 0.0_ccs_real
    self%payload%bP_f = 0.0_ccs_real
    self%payload%correction_f = 0.0_ccs_real
    self%payload%phi_values_f = 0.0_ccs_real
    self%payload%gradients_f = 0.0_ccs_real
    self%payload%rvecs_f = 0.0_ccs_real

  end subroutine ensure_scalar_transport_capacity

end submodule fv_equations_scalar_transport
