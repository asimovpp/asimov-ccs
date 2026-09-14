module fv_equations
#include "ccs_macros.inc"

  use constants, only: ndim
  use fv_kernels, only: advection_kernel, diffusion_kernel
  use kinds, only: ccs_int, ccs_real
  use types, only: cell_locator, field, matrix_values, vector_values

  implicit none

  private

  public :: equation
  public :: momentum_equation
  public :: poisson_equation
  public :: scalar_transport_equation

  !> Matrix row and RHS entry assembled by an equation apply step.
  type :: equation_row
    integer(ccs_int) :: global_row_index = 0_ccs_int
    integer(ccs_int) :: n_entries = 0_ccs_int
    integer(ccs_int), allocatable :: global_col_indices(:)
    real(ccs_real), allocatable :: coefficients(:)
    real(ccs_real) :: rhs = 0.0_ccs_real
    logical :: is_ready = .false.
  end type equation_row

  !> Row-local data required for one scalar transport equation row.
  type :: scalar_transport_payload
    integer(ccs_int) :: index_p = 0_ccs_int
    integer(ccs_int) :: global_index_p = 0_ccs_int
    integer(ccs_int) :: nnb = 0_ccs_int
    integer(ccs_int), allocatable :: global_indices_nb(:)
    logical, allocatable :: is_boundary_f(:)
    real(ccs_real), allocatable :: flux_f(:)
    real(ccs_real), allocatable :: diff_coeff_f(:)
    real(ccs_real), allocatable :: interpolation_factor_f(:)
    real(ccs_real), allocatable :: aPb_f(:)
    real(ccs_real), allocatable :: bP_f(:)
    real(ccs_real), allocatable :: correction_f(:)
    real(ccs_real), allocatable :: phi_values_f(:, :)
    real(ccs_real), allocatable :: gradients_f(:, :, :)
    real(ccs_real), allocatable :: rvecs_f(:, :, :)
  end type scalar_transport_payload

  !> Row-local momentum data: scalar transport coefficients plus pressure source inputs.
  type, extends(scalar_transport_payload) :: momentum_payload
    real(ccs_real) :: volume_p = 0.0_ccs_real
    real(ccs_real) :: pressure_gradient_p = 0.0_ccs_real
  end type momentum_payload

  !> Row-local data required for one pressure Poisson equation row.
  type :: pressure_poisson_payload
    integer(ccs_int) :: index_p = 0_ccs_int
    integer(ccs_int) :: global_index_p = 0_ccs_int
    integer(ccs_int) :: nnb = 0_ccs_int
    integer(ccs_int), allocatable :: global_indices_nb(:)
    logical, allocatable :: is_boundary_f(:)
    real(ccs_real), allocatable :: coeff_f(:)
    real(ccs_real), allocatable :: coeff_nb(:)
    real(ccs_real), allocatable :: rhs_f(:)
    real(ccs_real), allocatable :: aPb_f(:)
    real(ccs_real), allocatable :: bP_f(:)
  end type pressure_poisson_payload

  !> Abstract finite-volume equation interface.
  type, abstract :: equation
    type(equation_row), private :: row
  contains
    procedure(equation_apply_interface), deferred :: apply
    procedure :: prepare_row => equation_prepare_row
    procedure :: add_row_entry => equation_add_row_entry
    procedure :: set_rhs => equation_set_rhs
    procedure :: flush_row => equation_flush_row
    procedure :: flush_rhs => equation_flush_rhs
  end type equation

  abstract interface
    subroutine equation_apply_interface(self)
      import :: equation

      class(equation), intent(inout) :: self
    end subroutine equation_apply_interface
  end interface

  !> Scalar transport equation assembled from advection and diffusion kernels.
  type, extends(equation) :: scalar_transport_equation
    type(diffusion_kernel) :: diff_kernel
    class(advection_kernel), allocatable :: adv_kernel
    class(scalar_transport_payload), allocatable, private :: payload
    integer(ccs_int) :: capacity = 0_ccs_int
    integer(ccs_int) :: component = 0_ccs_int
  contains
    procedure :: init => scalar_transport_init
    procedure :: gather => scalar_transport_gather
    procedure :: apply => scalar_transport_apply
    procedure :: set_advection => scalar_transport_set_advection
  end type scalar_transport_equation

  !> Momentum equation; transport assembly is inherited from scalar transport.
  type, extends(scalar_transport_equation) :: momentum_equation
  contains
    procedure :: init => momentum_init
    procedure :: gather_pressure_source => momentum_gather_pressure_source
    procedure :: apply_pressure_source => momentum_apply_pressure_source
  end type momentum_equation

  !> Pressure Poisson equation assembled from diffusion kernels.
  type, extends(equation) :: poisson_equation
    type(diffusion_kernel) :: diff_kernel
    type(pressure_poisson_payload), private :: payload
    integer(ccs_int) :: capacity = 0_ccs_int
    logical :: fix_cached = .false.
    logical :: needs_fix = .false.
    integer(ccs_int) :: fix_row = -1_ccs_int
  contains
    procedure :: init => poisson_init
    procedure :: gather => poisson_gather
    procedure :: apply => poisson_apply
  end type poisson_equation

  interface
    module subroutine equation_prepare_row(self, global_row_index, required_capacity)
      class(equation), intent(inout) :: self
      integer(ccs_int), intent(in) :: global_row_index
      integer(ccs_int), intent(in), optional :: required_capacity
    end subroutine equation_prepare_row

    module subroutine equation_add_row_entry(self, global_col_index, coefficient)
      class(equation), intent(inout) :: self
      integer(ccs_int), intent(in) :: global_col_index
      real(ccs_real), intent(in) :: coefficient
    end subroutine equation_add_row_entry

    module subroutine equation_set_rhs(self, rhs)
      class(equation), intent(inout) :: self
      real(ccs_real), intent(in) :: rhs
    end subroutine equation_set_rhs

    module subroutine equation_flush_row(self, mat_coeffs, vec_values)
      class(equation), intent(in) :: self
      type(matrix_values), intent(inout) :: mat_coeffs
      type(vector_values), intent(inout) :: vec_values
    end subroutine equation_flush_row

    module subroutine equation_flush_rhs(self, vec_values)
      class(equation), intent(in) :: self
      type(vector_values), intent(inout) :: vec_values
    end subroutine equation_flush_rhs

    module subroutine ensure_row_capacity(row, required)
      type(equation_row), intent(inout) :: row
      integer(ccs_int), intent(in) :: required
    end subroutine ensure_row_capacity

    module subroutine scalar_transport_init(self, max_faces, component)
      class(scalar_transport_equation), intent(inout) :: self
      integer(ccs_int), intent(in) :: max_faces
      integer(ccs_int), intent(in), optional :: component
    end subroutine scalar_transport_init

    module subroutine scalar_transport_set_advection(self, kernel)
      class(scalar_transport_equation), intent(inout) :: self
      class(advection_kernel), intent(in) :: kernel
    end subroutine scalar_transport_set_advection

    module subroutine scalar_transport_gather(self, phi, loc_p, mass_flux, viscosity, density)
      class(scalar_transport_equation), intent(inout) :: self
      class(field), intent(in) :: phi
      type(cell_locator), intent(in) :: loc_p
      real(ccs_real), dimension(:), intent(in) :: mass_flux
      real(ccs_real), dimension(:), intent(in) :: viscosity
      real(ccs_real), dimension(:), intent(in) :: density
    end subroutine scalar_transport_gather

    module subroutine scalar_transport_apply(self)
      class(scalar_transport_equation), intent(inout) :: self
    end subroutine scalar_transport_apply

    module subroutine ensure_scalar_transport_capacity(self, required)
      class(scalar_transport_equation), intent(inout) :: self
      integer(ccs_int), intent(in) :: required
    end subroutine ensure_scalar_transport_capacity

    module subroutine momentum_init(self, max_faces, component)
      class(momentum_equation), intent(inout) :: self
      integer(ccs_int), intent(in) :: max_faces
      integer(ccs_int), intent(in), optional :: component
    end subroutine momentum_init

    module subroutine momentum_gather_pressure_source(self, p_gradient, loc_p)
      class(momentum_equation), intent(inout) :: self
      real(ccs_real), dimension(:), intent(in) :: p_gradient
      type(cell_locator), intent(in) :: loc_p
    end subroutine momentum_gather_pressure_source

    module subroutine momentum_apply_pressure_source(self)
      class(momentum_equation), intent(inout) :: self
    end subroutine momentum_apply_pressure_source

    module subroutine ensure_momentum_payload(self)
      class(momentum_equation), intent(inout) :: self
    end subroutine ensure_momentum_payload

    module subroutine poisson_init(self, max_faces)
      class(poisson_equation), intent(inout) :: self
      integer(ccs_int), intent(in) :: max_faces
    end subroutine poisson_init

    module subroutine poisson_gather(self, phi, loc_p, inv_diagonal)
      class(poisson_equation), intent(inout) :: self
      class(field), intent(in) :: phi
      type(cell_locator), intent(in) :: loc_p
      real(ccs_real), dimension(:), intent(in) :: inv_diagonal
    end subroutine poisson_gather

    module subroutine poisson_apply(self)
      class(poisson_equation), intent(inout) :: self
    end subroutine poisson_apply

    module subroutine ensure_poisson_fix(self, phi)
      class(poisson_equation), intent(inout) :: self
      class(field), intent(in) :: phi
    end subroutine ensure_poisson_fix

    module subroutine ensure_poisson_capacity(self, required)
      class(poisson_equation), intent(inout) :: self
      integer(ccs_int), intent(in) :: required
    end subroutine ensure_poisson_capacity
  end interface

end module fv_equations
