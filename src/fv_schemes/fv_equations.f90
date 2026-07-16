module fv_equations
#include "ccs_macros.inc"

  ! TODO: trim these
  use kinds, only: ccs_int, ccs_real
  use fv_kernels, only: diffusion_kernel, advection_kernel
  use types, only: vector_values, cell_locator, &
                   face_locator, neighbour_locator, matrix_values, field
  use constants, only: ndim
  use utils, only: set_entry, set_row, set_col, exit_print
  use bc_constants, only: bc_type_dirichlet
  use meshing, only: get_face_area, get_global_index, get_local_index, count_neighbours, &
                     get_boundary_status, get_face_normal, create_neighbour_locator, create_face_locator, &
                     create_cell_locator, get_volume, get_distance, &
                     get_face_interpolation, get_global_num_cells, get_centre

  implicit none

  public

  !v Abstract finite-volume equation interface
  !
  !  Coordinates the kernel-driven lifecycle described in theory/kernels.md by
  !  providing hooks for payload initialisation, data gathering, and assembly of
  !  matrix and RHS contributions at a cell.
  type, abstract :: equation
  contains
    procedure(equation_init_interface), deferred :: init
    procedure(equation_gather_interface), deferred :: gather
    procedure(equation_apply_interface), deferred :: apply
  end type equation

  abstract interface
    !> Prepare an equation-specific payload before sweeping over cells
    !> Switch to flow and pull out information
    !> Setup readonly views and add destructor
    subroutine equation_init_interface(self, max_faces, mf, viscosity, density, component)
      import :: equation
      import :: ccs_int
      import :: ccs_real

      class(equation), intent(inout) :: self
      integer(ccs_int), intent(in) :: max_faces                    !< Maximum stencil width exposed by neighbour faces
      real(ccs_real), dimension(:), target, intent(in), optional :: mf        !< Mass flux field used by advection kernels
      real(ccs_real), dimension(:), target, intent(in), optional :: viscosity !< Effective viscosity per control volume
      real(ccs_real), dimension(:), target, intent(in), optional :: density   !< Fluid density per control volume
      integer(ccs_int), intent(in), optional :: component                     !< Velocity component handled by the equation
    end subroutine equation_init_interface

    !> Gather the primitive data required by kernels for a given owner cell
    subroutine equation_gather_interface(self, phi, loc_p)
      import :: equation
      import :: field
      import :: cell_locator

      class(equation), intent(inout) :: self
      class(field), intent(in) :: phi       !< Field providing values and gradients
      type(cell_locator), intent(in) :: loc_p !< Locator for the owner cell being assembled
    end subroutine equation_gather_interface

    !> Apply kernel contributions for the current payload into the matrix row
    pure subroutine equation_apply_interface(self, mat_coeffs, vec_values)
      import :: equation
      import :: field
      import :: cell_locator
      import :: matrix_values
      import :: vector_values

      class(equation), intent(inout) :: self
      type(matrix_values), intent(inout) :: mat_coeffs !< Assembled implicit coefficients container
      type(vector_values), intent(inout) :: vec_values !< Assembled explicit/source contributions
    end subroutine equation_apply_interface
  end interface

  !v Face-oriented payload backing the pressure Poisson equation assembly
  type poisson_payload
    logical, allocatable :: is_boundary(:)
    integer(ccs_int) :: owner_local = 0
    integer(ccs_int) :: owner_global = 0
    integer(ccs_int) :: nfaces = 0
    integer(ccs_int), allocatable :: nb_global(:)
    real(ccs_real), allocatable :: coeff_face(:)
    real(ccs_real), allocatable :: coeff_nb(:)
    real(ccs_real), allocatable :: rhs_face(:)
    real(ccs_real), allocatable :: aPb(:)
    real(ccs_real), allocatable :: bP(:)
  end type poisson_payload

  !v Pressure Poisson equation constructed from diffusion kernels
  type, extends(equation) :: poisson_equation
    logical :: fix_cached = .false.
    logical :: needs_fix = .false.
    integer(ccs_int) :: fix_row = huge(0_ccs_int)
    integer(ccs_int) :: capacity = huge(0_ccs_int)
    real(ccs_real), pointer :: invA_data(:) => null()
    type(diffusion_kernel) :: diff_kernel
    type(poisson_payload) :: data
  contains
    procedure :: init => poisson_init
    procedure :: gather => poisson_gather
    procedure :: apply => poisson_apply
    procedure :: bind_inverse => poisson_bind_inverse
  end type poisson_equation

  !v Face-orientated payload backing the momentum equation gather/apply cycle
  !
  !  Stores neighbour connectivity, transport metrics, and intermediate kernel
  !  data so that `momentum_equation%apply` can evaluate diffusion and advection
  !  contributions without re-querying the mesh.
  type momentum_payload
    logical, allocatable :: is_boundary(:)
    integer(ccs_int) :: owner_local = huge(0_ccs_int)
    integer(ccs_int) :: owner_global = huge(0_ccs_int)
    integer(ccs_int) :: nfaces = huge(0_ccs_int)
    integer(ccs_int), allocatable :: nb_global(:)
    real(ccs_real), allocatable :: flux(:)
    real(ccs_real), allocatable :: diff_coeff(:)
    real(ccs_real), allocatable :: lf(:)
    real(ccs_real), allocatable :: aPb(:)
    real(ccs_real), allocatable :: bP(:)
    real(ccs_real), allocatable :: face_correction(:)
    real(ccs_real), allocatable :: phi_vals(:, :)
    real(ccs_real), allocatable :: grads(:, :, :)
    real(ccs_real), allocatable :: rvecs(:, :, :)
  end type momentum_payload

  !v Momentum equation assembled from diffusion and advection kernels
  !
  !  Implements the abstract equation interface by caching per-face payload
  !  data during `gather` and reusing it in `apply` to form the linearised
  !  momentum balance for a single velocity component.
  type, extends(equation) :: momentum_equation
    integer(ccs_int) :: component = huge(0_ccs_int)
    integer(ccs_int) :: capacity = huge(0_ccs_int)
    real(ccs_real), pointer :: mf(:) => null()
    real(ccs_real), pointer :: viscosity(:) => null()
    real(ccs_real), pointer :: density(:) => null()
    type(momentum_payload) :: data
    type(diffusion_kernel) :: diff_kernel
    class(advection_kernel), allocatable :: adv_kernel
  contains
    procedure :: init => momentum_init
    procedure :: gather => momentum_gather
    procedure :: apply => momentum_apply
    procedure :: set_advection => momentum_set_advection
  end type momentum_equation

contains

  subroutine poisson_init(self, max_faces, mf, viscosity, density, component)

    class(poisson_equation), intent(inout) :: self
    integer(ccs_int), intent(in) :: max_faces
    real(ccs_real), dimension(:), target, intent(in), optional :: mf
    real(ccs_real), dimension(:), target, intent(in), optional :: viscosity
    real(ccs_real), dimension(:), target, intent(in), optional :: density
    integer(ccs_int), intent(in), optional :: component

    associate (foo => mf); end associate
    associate (foo => viscosity); end associate
    associate (foo => density); end associate
    associate (foo => component); end associate

    call ensure_poisson_capacity(self, max_faces)

    self%data%nfaces = 0
    self%data%owner_local = 0
    self%data%owner_global = 0
    self%fix_cached = .false.
    self%needs_fix = .false.
    self%fix_row = -1
    self%capacity = 0

  end subroutine poisson_init

  subroutine poisson_bind_inverse(self, invA)

    class(poisson_equation), intent(inout) :: self
    real(ccs_real), dimension(:), target, intent(in) :: invA

    self%invA_data => invA

  end subroutine poisson_bind_inverse

  subroutine poisson_gather(self, phi, loc_p)

    use fv, only: compute_boundary_coeffs

    class(poisson_equation), intent(inout) :: self
    class(field), intent(in) :: phi
    type(cell_locator), intent(in) :: loc_p

    type(neighbour_locator) :: loc_nb
    type(face_locator) :: loc_f
    integer(ccs_int) :: index_p
    integer(ccs_int) :: global_index_p
    integer(ccs_int) :: nnb
    integer(ccs_int) :: j
    integer(ccs_int) :: index_nb
    integer(ccs_int) :: global_nb
    logical :: is_boundary
    type(cell_locator) :: loc_nb_cell
    real(ccs_real) :: Vp
    real(ccs_real) :: V_nb
    real(ccs_real) :: Vf
    real(ccs_real) :: invA_p
    real(ccs_real) :: invA_nb
    real(ccs_real) :: invA_f
    real(ccs_real) :: face_area
    real(ccs_real) :: coeff_f
    real(ccs_real), dimension(2) :: coeffs
    real(ccs_real) :: explicit_term
    real(ccs_real) :: aPb
    real(ccs_real) :: bP
    real(ccs_real) :: interpol
    real(ccs_real), dimension(ndim) :: face_normal
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
    real(ccs_real), dimension(3, 2) :: rvecs
    real(ccs_real), dimension(3, 2) :: grads
    real(ccs_real), dimension(2) :: phi_dummy

    call assert_poisson_inverse_bound(self)
    call ensure_poisson_fix(self, phi)

    call get_local_index(loc_p, index_p)
    call get_global_index(loc_p, global_index_p)
    call count_neighbours(loc_p, nnb)

    call ensure_poisson_capacity(self, nnb)

    self%data%owner_local = index_p
    self%data%owner_global = global_index_p
    self%data%nfaces = nnb

    call get_volume(loc_p, Vp)
    invA_p = self%invA_data(index_p)

    do j = 1, nnb
      call create_face_locator(index_p, j, loc_f)
      call get_face_area(loc_f, face_area)
      call get_boundary_status(loc_f, is_boundary)

      rvecs = 0.0_ccs_real
      grads = 0.0_ccs_real
      phi_dummy = 0.0_ccs_real
      phi_dummy(1) = phi%values_ro(index_p)

      if (.not. is_boundary) then
        call create_neighbour_locator(loc_p, j, loc_nb)
        call get_local_index(loc_nb, index_nb)
        call get_global_index(loc_nb, global_nb)
        self%data%nb_global(j) = global_nb

        call get_face_interpolation(loc_f, interpol)
        call get_distance(loc_p, loc_nb, dx)
        dxmag = sqrt(sum(dx**2))

        call get_volume(loc_nb, V_nb)
        Vf = interpol * Vp + (1.0_ccs_real - interpol) * V_nb

        if (phi%enable_cell_corrections) then
          call get_face_normal(loc_f, face_normal)
          call get_centre(loc_p, x_p)
          call create_cell_locator(index_nb, loc_nb_cell)
          call get_centre(loc_nb_cell, x_nb)
          call get_centre(loc_f, x_f)

          dx_orth = min(dot_product(x_f - x_p, face_normal), dot_product(x_nb - x_f, face_normal))
          dxmag = 2.0_ccs_real * dx_orth
          x_nb_prime = x_f + dx_orth * face_normal
          x_p_prime = x_f - dx_orth * face_normal

          grad_phi_p = [phi%x_gradients_ro(index_p), phi%y_gradients_ro(index_p), phi%z_gradients_ro(index_p)]
          grad_phi_nb = [phi%x_gradients_ro(index_nb), phi%y_gradients_ro(index_nb), phi%z_gradients_ro(index_nb)]

          grads(1:ndim, 1) = grad_phi_p
          grads(1:ndim, 2) = grad_phi_nb
          rvecs(1:ndim, 1) = x_p_prime - x_p
          rvecs(1:ndim, 2) = x_nb_prime - x_nb
        end if

        phi_dummy = [phi%values_ro(index_p), phi%values_ro(index_nb)]

        invA_nb = self%invA_data(index_nb)
        invA_f = interpol * invA_p + (1.0_ccs_real - interpol) * invA_nb
        aPb = 0.0_ccs_real
        bP = 0.0_ccs_real
      else
        self%data%nb_global(j) = -1
        interpol = 0.5_ccs_real
        call get_distance(loc_p, loc_f, dx)
        dxmag = 2.0_ccs_real * sqrt(sum(dx**2))
        Vf = Vp
        invA_f = invA_p
        call get_face_normal(loc_f, face_normal)
        call compute_boundary_coeffs(phi, 0, loc_p, loc_f, face_normal, aPb, bP)
        phi_dummy(2) = aPb * phi_dummy(1) + bP
      end if

      coeff_f = face_area * Vf * invA_f / dxmag
      coeffs = self%diff_kernel%eval_coeffs(coeff_f)
      explicit_term = self%diff_kernel%eval_explicit(phi_dummy, coeff_f, interpol, rvecs, grads)

      self%data%is_boundary(j) = is_boundary
      self%data%coeff_face(j) = coeffs(1)
      self%data%rhs_face(j) = explicit_term

      if (.not. is_boundary) then
        self%data%coeff_nb(j) = coeffs(2)
        self%data%aPb(j) = 0.0_ccs_real
        self%data%bP(j) = 0.0_ccs_real
      else
        self%data%coeff_nb(j) = coeffs(2)
        self%data%aPb(j) = aPb
        self%data%bP(j) = bP
      end if
    end do

  end subroutine poisson_gather

  pure subroutine poisson_apply(self, mat_coeffs, vec_values)

    class(poisson_equation), intent(inout) :: self
    type(matrix_values), intent(inout) :: mat_coeffs
    type(vector_values), intent(inout) :: vec_values

    real(ccs_real) :: diag
    real(ccs_real) :: rhs
    integer(ccs_int) :: j
    integer(ccs_int) :: row

    row = self%data%owner_global
    diag = 0.0_ccs_real
    rhs = 0.0_ccs_real

    call set_row(row, mat_coeffs)
    call set_row(row, vec_values)

    do j = 1, self%data%nfaces
      diag = diag + self%data%coeff_face(j)
      rhs = rhs + self%data%rhs_face(j)

      if (.not. self%data%is_boundary(j)) then
        call set_col(self%data%nb_global(j), mat_coeffs)
        call set_entry(self%data%coeff_nb(j), mat_coeffs)
      else
        diag = diag + self%data%aPb(j) * self%data%coeff_nb(j)
        rhs = rhs - self%data%coeff_nb(j) * self%data%bP(j)
      end if
    end do

    if (self%needs_fix .and. row == self%fix_row) then
      diag = diag + 1.0e30_ccs_real
    end if

    call set_col(row, mat_coeffs)
    call set_entry(diag, mat_coeffs)
    call set_entry(rhs, vec_values)

  end subroutine poisson_apply

  subroutine assert_poisson_inverse_bound(self)

    class(poisson_equation), intent(in) :: self

    if (.not. associated(self%invA_data)) then
      error stop "poisson_equation%gather: inverse diagonal not bound"
    end if

  end subroutine assert_poisson_inverse_bound

  subroutine ensure_poisson_fix(self, phi)

    class(poisson_equation), intent(inout) :: self
    class(field), intent(in) :: phi

    integer(ccs_int) :: global_num_cells
    integer(ccs_int) :: cps

    if (self%fix_cached) return

    self%needs_fix = .not. any(phi%bcs%bc_types(:) == bc_type_dirichlet)
    if (self%needs_fix) then
      call get_global_num_cells(global_num_cells)
      cps = int(sqrt(real(global_num_cells, kind=ccs_real)), ccs_int)
      self%fix_row = (cps / 2) * (1 + cps)
    else
      self%fix_row = -1
    end if
    self%fix_cached = .true.

  end subroutine ensure_poisson_fix

  subroutine ensure_poisson_capacity(self, required)

    class(poisson_equation), intent(inout) :: self
    integer(ccs_int), intent(in) :: required

    integer(ccs_int) :: allocate_size

    allocate_size = max(1_ccs_int, required)

    if (allocate_size <= self%capacity) return

    if (self%capacity > 0) then
      if (allocated(self%data%nb_global)) deallocate (self%data%nb_global)
      if (allocated(self%data%is_boundary)) deallocate (self%data%is_boundary)
      if (allocated(self%data%coeff_face)) deallocate (self%data%coeff_face)
      if (allocated(self%data%coeff_nb)) deallocate (self%data%coeff_nb)
      if (allocated(self%data%rhs_face)) deallocate (self%data%rhs_face)
      if (allocated(self%data%aPb)) deallocate (self%data%aPb)
      if (allocated(self%data%bP)) deallocate (self%data%bP)
    end if

    allocate (self%data%nb_global(allocate_size))
    allocate (self%data%is_boundary(allocate_size))
    allocate (self%data%coeff_face(allocate_size))
    allocate (self%data%coeff_nb(allocate_size))
    allocate (self%data%rhs_face(allocate_size))
    allocate (self%data%aPb(allocate_size))
    allocate (self%data%bP(allocate_size))

    self%capacity = allocate_size

    self%data%nb_global = 0_ccs_int
    self%data%is_boundary = .false.
    self%data%coeff_face = 0.0_ccs_real
    self%data%coeff_nb = 0.0_ccs_real
    self%data%rhs_face = 0.0_ccs_real
    self%data%aPb = 0.0_ccs_real
    self%data%bP = 0.0_ccs_real

  end subroutine ensure_poisson_capacity

  !> Bind transport properties and allocate payload storage for the momentum equation
  subroutine momentum_init(self, max_faces, mf, viscosity, density, component)

    class(momentum_equation), intent(inout) :: self
    integer(ccs_int), intent(in) :: max_faces                    !< Maximum face count encountered during assembly
    real(ccs_real), dimension(:), target, intent(in), optional :: mf        !< Face mass flux field supplied by the system
    real(ccs_real), dimension(:), target, intent(in), optional :: viscosity !< Effective viscosity field
    real(ccs_real), dimension(:), target, intent(in), optional :: density   !< Density field used for diffusion/Schmidt number
    integer(ccs_int), intent(in), optional :: component                     !< Velocity component owned by this equation

    if (.not. present(mf) .or. .not. present(viscosity) .or. .not. present(density) .or. .not. present(component)) then
      error stop "momentum_equation%init: missing required arguments"
    end if

    self%mf => mf
    self%viscosity => viscosity
    self%density => density
    if (present(component)) then
      self%component = component
    else
      self%component = 0
    end if
    self%capacity = 0

    call ensure_payload_capacity(self, max_faces)
    self%data%nfaces = 0
    self%data%owner_local = 0
    self%data%owner_global = 0

  end subroutine momentum_init

  !> Configure the advection kernel used when assembling the momentum equation
  subroutine momentum_set_advection(self, kernel)

    class(momentum_equation), intent(inout) :: self
    class(advection_kernel), intent(in) :: kernel !< Kernel instance providing advective coefficients

    if (allocated(self%adv_kernel)) then
      deallocate (self%adv_kernel)
    end if
    allocate (self%adv_kernel, source=kernel)

  end subroutine momentum_set_advection

  !> Collect the face-wise payload for a given owner cell prior to kernel evaluation
  subroutine momentum_gather(self, phi, loc_p)

    use fv, only: calc_diffusion_coeff, compute_boundary_coeffs

    class(momentum_equation), intent(inout) :: self
    class(field), intent(in) :: phi     !< Field providing values, gradients, and options for the sweep
    type(cell_locator), intent(in) :: loc_p !< Locator of the owner cell whose row is being assembled

    type(neighbour_locator) :: loc_nb
    type(face_locator) :: loc_f
    integer(ccs_int) :: index_p
    integer(ccs_int) :: global_index_p
    integer(ccs_int) :: nnb
    integer(ccs_int) :: j
    logical :: is_boundary
    integer(ccs_int) :: index_nb
    integer(ccs_int) :: global_nb
    type(cell_locator) :: loc_nb_cell
    integer(ccs_int) :: face_index
    real(ccs_real) :: phi_p
    real(ccs_real) :: phi_nb
    real(ccs_real) :: flux
    real(ccs_real) :: flux_mass
    real(ccs_real) :: face_area
    real(ccs_real) :: diff_coeff
    real(ccs_real) :: interpol
    real(ccs_real) :: aPb
    real(ccs_real) :: bP
    real(ccs_real) :: face_corr
    real(ccs_real), dimension(ndim) :: face_normal
    real(ccs_real), dimension(ndim) :: grad_phi_p
    real(ccs_real), dimension(ndim) :: grad_phi_nb
    real(ccs_real), dimension(ndim) :: x_p
    real(ccs_real), dimension(ndim) :: x_nb
    real(ccs_real), dimension(ndim) :: x_f
    real(ccs_real), dimension(ndim) :: x_nb_prime
    real(ccs_real), dimension(ndim) :: x_p_prime
    real(ccs_real) :: dx_orth
    real(ccs_real) :: SchmidtNo

    if (.not. allocated(self%adv_kernel)) then
      error stop "momentum_equation%gather: advection kernel not configured"
    end if

    call get_local_index(loc_p, index_p)
    call get_global_index(loc_p, global_index_p)
    call count_neighbours(loc_p, nnb)

    call ensure_payload_capacity(self, nnb)

    self%data%owner_local = index_p
    self%data%owner_global = global_index_p
    self%data%nfaces = nnb

    phi_p = phi%values_ro(index_p)
    SchmidtNo = phi%Schmidt

    do j = 1, nnb
      call create_neighbour_locator(loc_p, j, loc_nb)
      call get_boundary_status(loc_nb, is_boundary)
      call create_face_locator(index_p, j, loc_f)
      call get_local_index(loc_f, face_index)
      call get_face_area(loc_f, face_area)
      if (.not. is_boundary) then
        call get_face_interpolation(loc_f, interpol)
        self%data%lf(j) = interpol
      else
        self%data%lf(j) = 0.5_ccs_real
      end if


      self%data%is_boundary(j) = is_boundary
      self%data%grads(:, :, j) = 0.0_ccs_real
      self%data%rvecs(:, :, j) = 0.0_ccs_real
      self%data%aPb(j) = 0.0_ccs_real
      self%data%bP(j) = 0.0_ccs_real

      flux_mass = self%mf(face_index)

      face_corr = 0.0_ccs_real

      if (.not. is_boundary) then
        call get_local_index(loc_nb, index_nb)
        call get_global_index(loc_nb, global_nb)
        self%data%nb_global(j) = global_nb

        if (index_nb < index_p) flux_mass = -flux_mass

        phi_nb = phi%values_ro(index_nb)

        if (phi%enable_cell_corrections) then
          call get_face_normal(loc_f, face_normal)
          call get_centre(loc_p, x_p)
          call create_cell_locator(index_nb, loc_nb_cell)
          call get_centre(loc_nb_cell, x_nb)
          call get_centre(loc_f, x_f)

          dx_orth = min(dot_product(x_f - x_p, face_normal), dot_product(x_nb - x_f, face_normal))
          x_nb_prime = x_f + dx_orth * face_normal
          x_p_prime = x_f - dx_orth * face_normal

          grad_phi_p = [phi%x_gradients_ro(index_p), phi%y_gradients_ro(index_p), phi%z_gradients_ro(index_p)]
          grad_phi_nb = [phi%x_gradients_ro(index_nb), phi%y_gradients_ro(index_nb), phi%z_gradients_ro(index_nb)]

          self%data%grads(:, 1, j) = grad_phi_p
          self%data%grads(:, 2, j) = grad_phi_nb
          self%data%rvecs(:, 1, j) = x_p_prime - x_p
          self%data%rvecs(:, 2, j) = x_nb_prime - x_nb

          face_corr = 0.5_ccs_real * (dot_product(grad_phi_p, self%data%rvecs(:, 1, j)) + &
                                      dot_product(grad_phi_nb, self%data%rvecs(:, 2, j)))
        end if

        call calc_diffusion_coeff(index_p, j, phi%enable_cell_corrections, &
                                  self%viscosity(index_p), self%viscosity(index_nb), &
                                  self%density(index_p), self%density(index_nb), &
                                  SchmidtNo, diff_coeff)
      else
        self%data%nb_global(j) = -1_ccs_int

        call get_face_normal(loc_f, face_normal)
        call compute_boundary_coeffs(phi, self%component, loc_p, loc_f, face_normal, aPb, bP)

        self%data%aPb(j) = aPb
        self%data%bP(j) = bP

        phi_nb = aPb * phi_p + bP

        call calc_diffusion_coeff(index_p, j, .false., &
                                  self%viscosity(index_p), 0.0_ccs_real, &
                                  self%density(index_p), 0.0_ccs_real, &
                                  SchmidtNo, diff_coeff)
        diff_coeff = diff_coeff / 2.0_ccs_real
      end if

      flux = flux_mass * face_area
      self%data%flux(j) = flux
      self%data%diff_coeff(j) = diff_coeff
      self%data%phi_vals(1, j) = phi_p
      self%data%phi_vals(2, j) = phi_nb
      self%data%face_correction(j) = face_corr
    end do

  end subroutine momentum_gather

  !> Assemble momentum equation coefficients and RHS contributions from the cached payload
  pure subroutine momentum_apply(self, mat_coeffs, vec_values)

    class(momentum_equation), intent(inout) :: self
    type(matrix_values), intent(inout) :: mat_coeffs !< Matrix coefficient container to populate
    type(vector_values), intent(inout) :: vec_values !< RHS container to populate

    real(ccs_real), dimension(2) :: diff_coeffs
    real(ccs_real), dimension(2) :: adv_coeffs
    real(ccs_real) :: diff_expl
    real(ccs_real) :: adv_expl
    real(ccs_real) :: a_P
    real(ccs_real) :: a_F
    real(ccs_real) :: rhs
    integer(ccs_int) :: j
    integer(ccs_int) :: row

    row = self%data%owner_global
    a_P = 0.0_ccs_real
    rhs = 0.0_ccs_real

    call set_row(row, mat_coeffs)
    call set_row(row, vec_values)

    do j = 1, self%data%nfaces
      diff_coeffs = self%diff_kernel%eval_coeffs(self%data%diff_coeff(j))
      diff_expl = self%diff_kernel%eval_explicit(self%data%phi_vals(:, j), &
                                                     self%data%diff_coeff(j), &
                                                     self%data%lf(j), &
                                                     self%data%rvecs(:, :, j), &
                                                     self%data%grads(:, :, j))

      adv_coeffs = self%adv_kernel%eval_coeffs(self%data%flux(j))
      adv_expl = self%adv_kernel%eval_explicit(self%data%phi_vals(:, j), &
                                                   self%data%flux(j), &
                                                   self%data%lf(j), &
                                                   self%data%rvecs(:, :, j), &
                                                   self%data%grads(:, :, j))

      rhs = rhs - diff_expl - adv_expl - self%data%face_correction(j) * self%data%flux(j)

      ! Use the advection kernel's owner coefficient for the implicit contribution.
      a_P = a_P + adv_coeffs(1) + diff_coeffs(1)
      a_F = adv_coeffs(2) + diff_coeffs(2)

      if (.not. self%data%is_boundary(j)) then
        call set_col(self%data%nb_global(j), mat_coeffs)
        call set_entry(a_F, mat_coeffs)
      else
        a_P = a_P + self%data%aPb(j) * a_F
        rhs = rhs - a_F * self%data%bP(j)
      end if
    end do

    call set_col(row, mat_coeffs)
    call set_entry(a_P, mat_coeffs)
    call set_entry(rhs, vec_values)

  end subroutine momentum_apply

  !> Ensure the cached payload buffers can accommodate the requested number of faces
  subroutine ensure_payload_capacity(self, required)

    class(momentum_equation), intent(inout) :: self
    integer(ccs_int), intent(in) :: required !< Number of faces that must fit in the payload buffers
    integer(ccs_int) :: diff_width
    integer(ccs_int) :: adv_width
    integer(ccs_int) :: stencil_width

    diff_width = self%diff_kernel%get_width()
    if (allocated(self%adv_kernel)) then
      adv_width = self%adv_kernel%get_width()
    else
      adv_width = 1_ccs_int
    end if
    stencil_width = max(diff_width, adv_width) + 1_ccs_int

    if (stencil_width /= 2_ccs_int) then
      error stop "momentum_equation%ensure_payload_capacity: kernel width > 1 not yet supported"
    end if

    if (required <= self%capacity) return

    if (self%capacity > 0) then
      if (allocated(self%data%nb_global)) deallocate (self%data%nb_global)
      if (allocated(self%data%is_boundary)) deallocate (self%data%is_boundary)
      if (allocated(self%data%flux)) deallocate (self%data%flux)
      if (allocated(self%data%diff_coeff)) deallocate (self%data%diff_coeff)
      if (allocated(self%data%lf)) deallocate (self%data%lf)
      if (allocated(self%data%aPb)) deallocate (self%data%aPb)
      if (allocated(self%data%bP)) deallocate (self%data%bP)
      if (allocated(self%data%face_correction)) deallocate (self%data%face_correction)
      if (allocated(self%data%phi_vals)) deallocate (self%data%phi_vals)
      if (allocated(self%data%grads)) deallocate (self%data%grads)
      if (allocated(self%data%rvecs)) deallocate (self%data%rvecs)
    end if

    allocate (self%data%nb_global(required))
    allocate (self%data%is_boundary(required))
    allocate (self%data%flux(required))
    allocate (self%data%diff_coeff(required))
    allocate (self%data%lf(required))
    allocate (self%data%aPb(required))
    allocate (self%data%bP(required))
    allocate (self%data%face_correction(required))
    allocate (self%data%phi_vals(stencil_width, required))
    allocate (self%data%grads(ndim, stencil_width, required))
    allocate (self%data%rvecs(ndim, stencil_width, required))

    self%capacity = required

    self%data%nb_global = 0_ccs_int
    self%data%is_boundary = .false.
    self%data%flux = 0.0_ccs_real
    self%data%diff_coeff = 0.0_ccs_real
    self%data%lf = 0.0_ccs_real
    self%data%aPb = 0.0_ccs_real
    self%data%bP = 0.0_ccs_real
    self%data%face_correction = 0.0_ccs_real
    self%data%phi_vals = 0.0_ccs_real
    self%data%grads = 0.0_ccs_real
    self%data%rvecs = 0.0_ccs_real

  end subroutine ensure_payload_capacity

end module fv_equations
