!v Submodule file fv_equations_poisson.smod
submodule(fv_equations) fv_equations_poisson

  use bc_constants, only: bc_type_dirichlet
  use fv, only: compute_boundary_coeffs
  use meshing, only: count_neighbours, create_face_locator, create_neighbour_locator, &
                     get_boundary_status, get_centre, get_distance, get_face_area, get_face_interpolation, &
                     get_face_normal, get_global_index, get_global_num_cells, get_local_index, get_volume
  use types, only: face_locator, neighbour_locator

  implicit none

contains

  module subroutine poisson_init(self, max_faces)
    class(poisson_equation), intent(inout) :: self
    integer(ccs_int), intent(in) :: max_faces

    call ensure_poisson_capacity(self, max_faces)
    call ensure_row_capacity(self%row, max_faces + 1_ccs_int)

    self%payload%nnb = 0_ccs_int
    self%payload%index_p = 0_ccs_int
    self%payload%global_index_p = 0_ccs_int
    self%fix_cached = .false.
    self%needs_fix = .false.
    self%fix_row = -1_ccs_int

  end subroutine poisson_init

  module subroutine poisson_gather(self, phi, loc_p, inv_diagonal)
    class(poisson_equation), intent(inout) :: self
    class(field), intent(in) :: phi
    type(cell_locator), intent(in) :: loc_p
    real(ccs_real), dimension(:), intent(in) :: inv_diagonal

    type(neighbour_locator) :: loc_nb
    type(face_locator) :: loc_f
    integer(ccs_int) :: index_p
    integer(ccs_int) :: global_index_p
    integer(ccs_int) :: nnb
    integer(ccs_int) :: j
    integer(ccs_int) :: index_nb
    integer(ccs_int) :: global_index_nb
    logical :: is_boundary_f
    real(ccs_real) :: volume_p
    real(ccs_real) :: volume_nb
    real(ccs_real) :: volume_f
    real(ccs_real) :: inv_diagonal_p
    real(ccs_real) :: inv_diagonal_nb
    real(ccs_real) :: inv_diagonal_f
    real(ccs_real) :: area_f
    real(ccs_real) :: coeff_f
    real(ccs_real), dimension(2) :: coeffs
    real(ccs_real) :: explicit_term_f
    real(ccs_real) :: aPb_f
    real(ccs_real) :: bP_f
    real(ccs_real) :: interpolation_factor_f
    real(ccs_real), dimension(ndim) :: normal_f
    real(ccs_real), dimension(ndim) :: dx
    real(ccs_real) :: dxmag
    real(ccs_real) :: dx_orth
    real(ccs_real), dimension(ndim) :: x_p
    real(ccs_real), dimension(ndim) :: x_nb
    real(ccs_real), dimension(ndim) :: x_f
    real(ccs_real), dimension(ndim) :: x_nb_prime
    real(ccs_real), dimension(ndim) :: x_p_prime
    real(ccs_real), dimension(ndim) :: grad_phi_p
    real(ccs_real), dimension(ndim) :: grad_phi_nb
    real(ccs_real), dimension(3, 2) :: rvecs_f
    real(ccs_real), dimension(3, 2) :: gradients_f
    real(ccs_real), dimension(2) :: phi_values_f

    call ensure_poisson_fix(self, phi)

    call get_local_index(loc_p, index_p)
    call get_global_index(loc_p, global_index_p)
    call count_neighbours(loc_p, nnb)

    call ensure_poisson_capacity(self, nnb)

    self%payload%index_p = index_p
    self%payload%global_index_p = global_index_p
    self%payload%nnb = nnb

    call get_volume(loc_p, volume_p)
    inv_diagonal_p = inv_diagonal(index_p)

    do j = 1, nnb
      call create_face_locator(index_p, j, loc_f)
      call get_face_area(loc_f, area_f)
      call get_boundary_status(loc_f, is_boundary_f)

      rvecs_f = 0.0_ccs_real
      gradients_f = 0.0_ccs_real
      phi_values_f = 0.0_ccs_real
      phi_values_f(1) = phi%values_ro(index_p)

      if (.not. is_boundary_f) then
        call create_neighbour_locator(loc_p, j, loc_nb)
        call get_local_index(loc_nb, index_nb)
        call get_global_index(loc_nb, global_index_nb)
        self%payload%global_indices_nb(j) = global_index_nb

        call get_face_interpolation(loc_f, interpolation_factor_f)
        call get_distance(loc_p, loc_nb, dx)
        dxmag = sqrt(sum(dx**2))

        call get_volume(loc_nb, volume_nb)
        volume_f = interpolation_factor_f * volume_p + (1.0_ccs_real - interpolation_factor_f) * volume_nb

        if (phi%enable_cell_corrections) then
          call get_face_normal(loc_f, normal_f)
          call get_centre(loc_p, x_p)
          call get_centre(loc_nb, x_nb)
          call get_centre(loc_f, x_f)

          dx_orth = min(dot_product(x_f - x_p, normal_f), dot_product(x_nb - x_f, normal_f))
          dxmag = 2.0_ccs_real * dx_orth
          x_nb_prime = x_f + dx_orth * normal_f
          x_p_prime = x_f - dx_orth * normal_f

          grad_phi_p = [phi%x_gradients_ro(index_p), phi%y_gradients_ro(index_p), phi%z_gradients_ro(index_p)]
          grad_phi_nb = [phi%x_gradients_ro(index_nb), phi%y_gradients_ro(index_nb), phi%z_gradients_ro(index_nb)]

          gradients_f(1:ndim, 1) = grad_phi_p
          gradients_f(1:ndim, 2) = grad_phi_nb
          rvecs_f(1:ndim, 1) = x_p_prime - x_p
          rvecs_f(1:ndim, 2) = x_nb_prime - x_nb
        end if

        phi_values_f = [phi%values_ro(index_p), phi%values_ro(index_nb)]

        inv_diagonal_nb = inv_diagonal(index_nb)
        inv_diagonal_f = interpolation_factor_f * inv_diagonal_p + &
                         (1.0_ccs_real - interpolation_factor_f) * inv_diagonal_nb
        aPb_f = 0.0_ccs_real
        bP_f = 0.0_ccs_real
      else
        self%payload%global_indices_nb(j) = -1_ccs_int
        interpolation_factor_f = 0.5_ccs_real
        call get_distance(loc_p, loc_f, dx)
        dxmag = 2.0_ccs_real * sqrt(sum(dx**2))
        volume_f = volume_p
        inv_diagonal_f = inv_diagonal_p
        call get_face_normal(loc_f, normal_f)
        call compute_boundary_coeffs(phi, 0_ccs_int, loc_p, loc_f, normal_f, aPb_f, bP_f)
        phi_values_f(2) = aPb_f * phi_values_f(1) + bP_f
      end if

      coeff_f = area_f * volume_f * inv_diagonal_f / dxmag
      coeffs = self%diff_kernel%eval_coeffs(coeff_f)
      explicit_term_f = self%diff_kernel%eval_explicit(phi_values_f, coeff_f, interpolation_factor_f, &
                                                       rvecs_f, gradients_f)

      self%payload%is_boundary_f(j) = is_boundary_f
      self%payload%coeff_f(j) = coeffs(1)
      self%payload%rhs_f(j) = explicit_term_f

      if (.not. is_boundary_f) then
        self%payload%coeff_nb(j) = coeffs(2)
        self%payload%aPb_f(j) = 0.0_ccs_real
        self%payload%bP_f(j) = 0.0_ccs_real
      else
        self%payload%coeff_nb(j) = coeffs(2)
        self%payload%aPb_f(j) = aPb_f
        self%payload%bP_f(j) = bP_f
      end if
    end do

  end subroutine poisson_gather

  module subroutine poisson_apply(self)
    class(poisson_equation), intent(inout) :: self

    real(ccs_real) :: diag
    real(ccs_real) :: rhs
    integer(ccs_int) :: j
    integer(ccs_int) :: global_row_index

    global_row_index = self%payload%global_index_p
    diag = 0.0_ccs_real
    rhs = 0.0_ccs_real

    call self%prepare_row(global_row_index, self%payload%nnb + 1_ccs_int)

    do j = 1, self%payload%nnb
      diag = diag + self%payload%coeff_f(j)
      rhs = rhs + self%payload%rhs_f(j)

      if (.not. self%payload%is_boundary_f(j)) then
        call self%add_row_entry(self%payload%global_indices_nb(j), self%payload%coeff_nb(j))
      else
        diag = diag + self%payload%aPb_f(j) * self%payload%coeff_nb(j)
        rhs = rhs - self%payload%coeff_nb(j) * self%payload%bP_f(j)
      end if
    end do

    if (self%needs_fix .and. global_row_index == self%fix_row) then
      diag = diag + 1.0e30_ccs_real
    end if

    call self%add_row_entry(global_row_index, diag)
    call self%set_rhs(rhs)

  end subroutine poisson_apply

  module subroutine ensure_poisson_fix(self, phi)
    class(poisson_equation), intent(inout) :: self
    class(field), intent(in) :: phi

    integer(ccs_int) :: global_num_cells
    integer(ccs_int) :: cps

    if (self%fix_cached) return

    self%needs_fix = .not. any(phi%bcs%bc_types(:) == bc_type_dirichlet)
    if (self%needs_fix) then
      call get_global_num_cells(global_num_cells)
      cps = int(sqrt(real(global_num_cells, ccs_real)), ccs_int)
      self%fix_row = (cps / 2) * (1 + cps)
    else
      self%fix_row = -1_ccs_int
    end if
    self%fix_cached = .true.

  end subroutine ensure_poisson_fix

  module subroutine ensure_poisson_capacity(self, required)
    class(poisson_equation), intent(inout) :: self
    integer(ccs_int), intent(in) :: required

    integer(ccs_int) :: allocate_size

    allocate_size = max(1_ccs_int, required)
    if (allocate_size <= self%capacity) return

    if (allocated(self%payload%global_indices_nb)) deallocate (self%payload%global_indices_nb)
    if (allocated(self%payload%is_boundary_f)) deallocate (self%payload%is_boundary_f)
    if (allocated(self%payload%coeff_f)) deallocate (self%payload%coeff_f)
    if (allocated(self%payload%coeff_nb)) deallocate (self%payload%coeff_nb)
    if (allocated(self%payload%rhs_f)) deallocate (self%payload%rhs_f)
    if (allocated(self%payload%aPb_f)) deallocate (self%payload%aPb_f)
    if (allocated(self%payload%bP_f)) deallocate (self%payload%bP_f)

    allocate (self%payload%global_indices_nb(allocate_size))
    allocate (self%payload%is_boundary_f(allocate_size))
    allocate (self%payload%coeff_f(allocate_size))
    allocate (self%payload%coeff_nb(allocate_size))
    allocate (self%payload%rhs_f(allocate_size))
    allocate (self%payload%aPb_f(allocate_size))
    allocate (self%payload%bP_f(allocate_size))

    self%capacity = allocate_size

    self%payload%global_indices_nb = 0_ccs_int
    self%payload%is_boundary_f = .false.
    self%payload%coeff_f = 0.0_ccs_real
    self%payload%coeff_nb = 0.0_ccs_real
    self%payload%rhs_f = 0.0_ccs_real
    self%payload%aPb_f = 0.0_ccs_real
    self%payload%bP_f = 0.0_ccs_real

  end subroutine ensure_poisson_capacity

end submodule fv_equations_poisson
