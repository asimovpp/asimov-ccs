!> Verify that optional residual storage preserves coupled and scalar-only solutions.
!>
!> Storage-off and storage-on cases run in separate processes. Each case advances two timesteps on a four-cell-per-side
!> cube. The coupled comparison covers all three velocity components, pressure and its
!> correction, scalar concentration, corrected cell mass imbalance and residual history.
!> The scalar-only case exercises the solver path without an optional residual monitor
!> and checks that each timestep replaces every stored scalar residual sentinel.

program test_residual_storage
#include "ccs_macros.inc"
  use core
  use kinds, only: ccs_int, ccs_real, CCS_MPI_PRECISION
  use types, only: fluid, field, cell_locator, face_locator, neighbour_locator, ccs_vector
  use parallel, only: initialise_parallel_environment, cleanup_parallel_environment, is_root
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi
  use fields, only: get_field, dealloc_fluid_fields
  use meshing, only: create_cell_locator, create_face_locator, create_neighbour_locator, &
                     get_centre, get_natural_index, get_local_num_cells, count_neighbours, &
                     get_local_index, get_face_area, get_boundary_status, nullify_mesh_object
  use vec, only: get_vector_data, restore_vector_data, get_vector_data_readonly, restore_vector_data_readonly
  use profiler, only: profiler_init, profiler_shutdown, profiler_begin_region, profiler_end_region
  use testing_lib, only: assert_bool, assert_eq, assert_gt, assert_close, par_env
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: iostat_end
  use mpi

  implicit none

  integer(ccs_int), parameter :: cells_per_side = 4, num_cells = cells_per_side**3
  integer(ccs_int), parameter :: num_steps = 2, num_iters = 3
  integer(ccs_int), parameter :: u_column = 1, v_column = 2, w_column = 3
  integer(ccs_int), parameter :: p_column = 4, p_prime_column = 5, c_column = 6
  character(len=7), dimension(*), parameter :: field_names = [character(len=7) :: "u", "v", "w", "p", "p_prime", "c"]
  integer(ccs_int), parameter :: imbalance_column = size(field_names) + 1
  integer(ccs_int), parameter :: pressure_residual_column = 1, v_change_column = 2, w_change_column = 3
  integer(ccs_int), parameter :: num_diagnostics = 3
  integer(ccs_int), parameter :: num_history_norms = 2 * (size(field_names) - 1)
  real(ccs_real), parameter :: scalar_residual_sentinel = 1.23456789e30_ccs_real
  real(ccs_real), parameter :: state_rtol = 1.0e-10_ccs_real, atol = 1.0e-12_ccs_real
  ! The production e12.4 history format rounds norms to four significant digits.
  real(ccs_real), parameter :: history_rtol = 2.0e-4_ccs_real
  real(ccs_real), parameter :: u_amplitude = 0.1_ccs_real, v_amplitude = 0.2_ccs_real, w_amplitude = 0.3_ccs_real
  class(parallel_environment), allocatable, target :: shared_env
  type(ccs_options) :: run_options
  type(fluid) :: flow
  class(field), pointer :: phi
  character(len=1024) :: mode, off_path, on_path
  integer :: state_unit
  integer(ccs_int) :: field_column, completed_step
  logical :: scalar_only

  call initialise_parallel_environment(par_env)
  call get_command_argument(1, mode)
  if (trim(mode) == "compare" .or. trim(mode) == "compare_scalar") then
    call get_command_argument(2, off_path)
    call get_command_argument(3, on_path)
    if (is_root(par_env)) then
      call compare_states(trim(off_path), trim(on_path))
      if (trim(mode) == "compare") call compare_history(trim(off_path), trim(on_path))
    end if
  else
    call get_environment_variable("CCS_RESIDUAL_SCALAR", mode)
    scalar_only = trim(mode) == "1"
    call profiler_init()
    call get_config(par_env, run_options)
    call configure_parallelism(run_options, par_env, shared_env)
    call profiler_begin_region("Total elapsed time")
    call profiler_begin_region("Total initialisation")
    call initialise_mesh(par_env, shared_env, run_options)
    call initialise_fields(par_env, run_options, flow)
    call initialise_flow(par_env, run_options, flow, get_init_flow, get_init_mass_flux)
    if (scalar_only) then
      do field_column = u_column, p_prime_column
        call get_field(flow, trim(field_names(field_column)), phi)
        phi%solver_parameters%solve = .false.
      end do
      call initialise_scalar_residual_sentinel(flow)
    end if
    completed_step = 0
    if (is_root(par_env)) open (newunit=state_unit, file="state.dat", status="replace")
    call profiler_end_region("Total initialisation")
    call run_solver(par_env, run_options, eval_sources, record_completed_step, flow)
    if (is_root(par_env)) close (state_unit)
    call profiler_end_region("Total elapsed time")
    call profiler_shutdown(par_env)
    call dealloc_fluid_fields(flow)
    call nullify_mesh_object()
  end if
  call cleanup_parallel_environment(par_env)

contains

  !> Initialise distinct velocity components with nonzero divergence on [0, pi] cubed.
  !>
  !> Different amplitudes expose component-specific errors, while the divergent initial
  !> velocity requires a nonzero pressure correction during the coupled solve.
  pure subroutine get_init_flow(loc_p, field_name, init_val)
    type(cell_locator), intent(in) :: loc_p !< Cell at which the analytic initial value is evaluated.
    character(len=*), intent(in) :: field_name !< Name selecting velocity, pressure or concentration.
    real(ccs_real), intent(inout) :: init_val !< Initial field value supplied to the solver callback.
    real(ccs_real), dimension(3) :: x

    call get_centre(loc_p, x)
    select case (field_name)
    case ("u")
      init_val = u_amplitude * sin(2 * x(1)) * sin(x(2)) * sin(x(3))
    case ("v")
      init_val = v_amplitude * sin(x(1)) * sin(2 * x(2)) * sin(x(3))
    case ("w")
      init_val = w_amplitude * sin(x(1)) * sin(x(2)) * sin(2 * x(3))
    case ("c")
      init_val = 1.0_ccs_real + 0.1_ccs_real * product(sin(x))
    case ("p", "p_prime")
      init_val = 0.0_ccs_real
    end select
  end subroutine get_init_flow

  !> Retain the core interpolation of the three initial velocity components onto faces.
  pure subroutine get_init_mass_flux(loc_f, init_val)
    type(face_locator), intent(in) :: loc_f !< Face supplied by the initialisation callback interface.
    real(ccs_real), intent(inout) :: init_val !< Interpolated mass flux retained without modification.

    associate (unused => loc_f, unused_value => init_val)
    end associate
  end subroutine get_init_mass_flux

  !> Supply zero integrated sources through the standard case callback.
  subroutine eval_sources(flow, phi, R, S)
    use fv, only: zero_sources
    type(fluid), intent(in) :: flow !< Fields used to evaluate the source terms.
    class(field), intent(in) :: phi !< Transported field whose sources are requested.
    class(ccs_vector), intent(inout) :: R !< Explicit integrated source, reset to zero.
    class(ccs_vector), intent(inout) :: S !< Implicit integrated source, reset to zero.

    call zero_sources(flow, phi, R, S)
  end subroutine eval_sources

  !> Fill allocated scalar residual storage with a large finite sentinel.
  !>
  !> The sentinel is reset before every scalar-only timestep, including the first.
  !> Checking its replacement after each step detects stale residuals even when the
  !> optional residual monitor is absent. A finite value permits an explicit magnitude
  !> check in addition to the independent finiteness assertion.
  subroutine initialise_scalar_residual_sentinel(flow)
    type(fluid), intent(in) :: flow !< Scalar field whose optional residual storage is initialised.
    class(field), pointer :: scalar
    real(ccs_real), dimension(:), pointer :: residual_values

    call get_field(flow, trim(field_names(c_column)), scalar)
    if (allocated(scalar%residuals)) then
      call get_vector_data(scalar%residuals, residual_values)
      residual_values = scalar_residual_sentinel
      call restore_vector_data(scalar%residuals, residual_values)
    end if
  end subroutine initialise_scalar_residual_sentinel

  !> Assemble and record one completed timestep in natural cell order.
  !>
  !> Owned cells contribute once to the state sum; diagnostics use global maxima.
  !> The allocatable environment is required by the core callback.
  subroutine record_completed_step(par_env, flow_fields)
    class(parallel_environment), allocatable, intent(in) :: par_env !< Solver communicator used for global assembly.
    type(fluid), intent(in) :: flow_fields !< Completed fields and optional stored residuals.
    real(ccs_real), dimension(size(field_names) + 1, num_cells) :: state_values
    real(ccs_real), dimension(num_diagnostics) :: diagnostics
    integer(ccs_int) :: natural_index_p
    integer :: ierr

    completed_step = completed_step + 1
    call collect_field_state(flow_fields, state_values, diagnostics)
    call collect_corrected_mass_imbalance(flow_fields, state_values(imbalance_column, :))
    call assert_scalar_residuals(flow_fields)
    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_Allreduce(MPI_IN_PLACE, state_values, size(state_values), CCS_MPI_PRECISION, MPI_SUM, par_env%comm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, diagnostics, size(diagnostics), CCS_MPI_PRECISION, MPI_MAX, par_env%comm, ierr)
    end select
    call assert_completed_state(flow_fields, state_values, diagnostics)
    if (is_root(par_env)) then
      do natural_index_p = 1, num_cells
        write (state_unit, '(2(i0,1x),*(es25.17e3,1x))') completed_step, natural_index_p, state_values(:, natural_index_p)
      end do
    end if
    if (scalar_only) call initialise_scalar_residual_sentinel(flow_fields)
  end subroutine record_completed_step

  !> Collect fields in natural cell order and maxima over owned cells.
  subroutine collect_field_state(flow_fields, state_values, diagnostics)
    type(fluid), intent(in) :: flow_fields !< Completed fields whose owned values are collected.
    real(ccs_real), dimension(:, :), intent(out) :: state_values !< Fields and a reserved corrected-imbalance column.
    real(ccs_real), dimension(:), intent(out) :: diagnostics !< Local maxima of pressure residual magnitude and changes in v and w.
    real(ccs_real), dimension(:), pointer :: field_values
    real(ccs_real) :: initial_value
    class(field), pointer :: current
    type(cell_locator) :: loc_p
    integer(ccs_int) :: nlocal, index_p, natural_index_p, column, diagnostic_column

    state_values = 0.0_ccs_real
    diagnostics = 0.0_ccs_real
    call get_local_num_cells(nlocal)
    do column = 1, size(field_names)
      call get_field(flow_fields, trim(field_names(column)), current)
      call get_vector_data_readonly(current%values, field_values)
      do index_p = 1, nlocal
        call create_cell_locator(index_p, loc_p)
        call get_natural_index(loc_p, natural_index_p)
        state_values(column, natural_index_p) = field_values(index_p)
        if (column == v_column .or. column == w_column) then
          initial_value = 0.0_ccs_real
          call get_init_flow(loc_p, trim(field_names(column)), initial_value)
          if (column == v_column) diagnostic_column = v_change_column
          if (column == w_column) diagnostic_column = w_change_column
          diagnostics(diagnostic_column) = max(diagnostics(diagnostic_column), abs(field_values(index_p) - initial_value))
        end if
      end do
      call restore_vector_data_readonly(current%values, field_values)
    end do
    call get_field(flow_fields, trim(field_names(p_column)), current)
    if (allocated(current%residuals) .and. .not. scalar_only) then
      call get_vector_data_readonly(current%residuals, field_values)
      diagnostics(pressure_residual_column) = max(0.0_ccs_real, maxval(abs(field_values(1:nlocal))))
      call restore_vector_data_readonly(current%residuals, field_values)
    end if
  end subroutine collect_field_state

  !> Integrate the oriented, face mass flux for every owned cell.
  subroutine collect_corrected_mass_imbalance(flow_fields, cell_imbalance)
    type(fluid), intent(in) :: flow_fields !< Completed flow containing the corrected face mass flux.
    real(ccs_real), dimension(:), intent(out) :: cell_imbalance !< Imbalance in natural cell order; zero for cells owned elsewhere.
    real(ccs_real), dimension(:), pointer :: face_flux
    real(ccs_real) :: oriented_area
    class(field), pointer :: mass_flux
    type(cell_locator) :: loc_p
    type(face_locator) :: loc_f
    type(neighbour_locator) :: loc_nb
    integer(ccs_int) :: nlocal, index_p, natural_index_p, index_f, index_nb, nnb, neighbour
    logical :: boundary

    cell_imbalance = 0.0_ccs_real
    call get_local_num_cells(nlocal)
    call get_field(flow_fields, "mf", mass_flux)
    call get_vector_data_readonly(mass_flux%values, face_flux)
    do index_p = 1, nlocal
      call create_cell_locator(index_p, loc_p)
      call get_natural_index(loc_p, natural_index_p)
      call count_neighbours(loc_p, nnb)
      do neighbour = 1, nnb
        call create_face_locator(index_p, neighbour, loc_f)
        call get_local_index(loc_f, index_f)
        call get_face_area(loc_f, oriented_area)
        call get_boundary_status(loc_f, boundary)
        if (.not. boundary) then
          call create_neighbour_locator(loc_p, neighbour, loc_nb)
          call get_local_index(loc_nb, index_nb)
          if (index_nb < index_p) oriented_area = -oriented_area
        end if
        cell_imbalance(natural_index_p) = cell_imbalance(natural_index_p) + oriented_area * face_flux(index_f)
      end do
    end do
    call restore_vector_data_readonly(mass_flux%values, face_flux)
  end subroutine collect_corrected_mass_imbalance

  !> Check that every owned stored scalar residual is finite and replaces the sentinel.
  subroutine assert_scalar_residuals(flow_fields)
    type(fluid), intent(in) :: flow_fields !< Completed scalar field with optional residual storage.
    class(field), pointer :: scalar
    real(ccs_real), dimension(:), pointer :: residual_values
    integer(ccs_int) :: nlocal

    call get_local_num_cells(nlocal)
    call get_field(flow_fields, trim(field_names(c_column)), scalar)
    if (allocated(scalar%residuals)) then
      call get_vector_data_readonly(scalar%residuals, residual_values)
      call assert_bool(all(ieee_is_finite(residual_values(1:nlocal))), "Stored scalar residuals must be finite")
      call assert_bool(all(abs(residual_values(1:nlocal)) < scalar_residual_sentinel / 2), &
                       "Scalar solve must overwrite every sentinel")
      call restore_vector_data_readonly(scalar%residuals, residual_values)
    end if
  end subroutine assert_scalar_residuals

  !> Verify finite state and a nontrivial coupled pressure and velocity update.
  subroutine assert_completed_state(flow_fields, state_values, diagnostics)
    type(fluid), intent(in) :: flow_fields !< Fields used to determine whether pressure residual storage exists.
    real(ccs_real), dimension(:, :), intent(in) :: state_values !< Globally assembled fields and corrected mass imbalance.
    real(ccs_real), dimension(:), intent(in) :: diagnostics !< Global maxima of pressure residual magnitude and velocity changes.
    class(field), pointer :: pressure

    call assert_bool(all(ieee_is_finite(state_values)), "All state values must be finite")
    if (.not. scalar_only) then
      call assert_gt(maxval(abs(state_values(p_prime_column, :))), atol, "Pressure correction must be nonzero")
      call assert_gt(diagnostics(v_change_column), atol, "v must evolve")
      call assert_gt(diagnostics(w_change_column), atol, "w must evolve")
      call get_field(flow_fields, trim(field_names(p_column)), pressure)
      if (allocated(pressure%residuals)) then
        call assert_gt(diagnostics(pressure_residual_column), atol, "Stored pressure residual must be nonzero")
      end if
    end if
  end subroutine assert_completed_state

  !> Compare every completed step and natural cell between the two independent runs.
  subroutine compare_states(off_path, on_path)
    character(len=*), intent(in) :: off_path !< Directory containing the storage-off state file.
    character(len=*), intent(in) :: on_path !< Directory containing the storage-on state file.
    integer :: off_unit, on_unit, ios
    integer(ccs_int) :: expected_step, expected_cell, off_step, on_step, off_cell, on_cell
    real(ccs_real), dimension(size(field_names) + 1) :: off_state, on_state

    open (newunit=off_unit, file=off_path // "/state.dat", status="old")
    open (newunit=on_unit, file=on_path // "/state.dat", status="old")
    do expected_step = 1, num_steps
      do expected_cell = 1, num_cells
        read (off_unit, *) off_step, off_cell, off_state
        read (on_unit, *) on_step, on_cell, on_state
        call assert_eq(off_step, expected_step, "Expected snapshot step")
        call assert_eq(off_cell, expected_cell, "Natural cell order")
        call assert_eq(on_step, off_step, "Same steps")
        call assert_eq(on_cell, off_cell, "Same natural cells")
        call assert_close(on_state, off_state, state_rtol, atol, &
                          "Storage must preserve fields and corrected mass imbalance")
      end do
    end do
    read (off_unit, *, iostat=ios) off_state
    call assert_eq(ios, iostat_end, "No additional storage-off snapshots")
    read (on_unit, *, iostat=ios) on_state
    call assert_eq(ios, iostat_end, "No additional storage-on snapshots")
    close (off_unit)
    close (on_unit)
  end subroutine compare_states

  !> Compare typed step, time, iteration and ten norms for each coupled iteration.
  !>
  !> History labels the step being advanced from zero; state output labels completed
  !> steps from one. The norm tolerance accommodates the production e12.4 formatting.
  subroutine compare_history(off_path, on_path)
    character(len=*), intent(in) :: off_path !< Directory containing storage-off residual history.
    character(len=*), intent(in) :: on_path !< Directory containing storage-on residual history.
    integer :: off_unit, on_unit, ios
    integer(ccs_int) :: expected_step, expected_iteration, off_step, on_step, off_iteration, on_iteration
    real(ccs_real) :: off_time, on_time
    real(ccs_real), dimension(num_history_norms) :: off_norms, on_norms

    open (newunit=off_unit, file=off_path // "/residuals.log", status="old")
    open (newunit=on_unit, file=on_path // "/residuals.log", status="old")
    ! Skip the residual-history headers.
    read (off_unit, '(a)')
    read (on_unit, '(a)')
    do expected_step = 0, num_steps - 1
      do expected_iteration = 1, num_iters
        read (off_unit, *) off_step, off_time, off_iteration, off_norms
        read (on_unit, *) on_step, on_time, on_iteration, on_norms
        call assert_eq(off_step, expected_step, "Expected history step")
        call assert_eq(off_iteration, expected_iteration, "Expected history iteration")
        call assert_eq(on_step, off_step, "Same history step")
        call assert_eq(on_time, off_time, "Same history time")
        call assert_eq(on_iteration, off_iteration, "Same history iteration")
        call assert_bool(all(ieee_is_finite(off_norms)), "Storage-off residual norms must be finite")
        call assert_bool(all(ieee_is_finite(on_norms)), "Storage-on residual norms must be finite")
        call assert_close(on_norms, off_norms, history_rtol, atol, "Storage must preserve logged norms")
      end do
    end do
    read (off_unit, *, iostat=ios) off_step, off_time, off_iteration, off_norms
    call assert_eq(ios, iostat_end, "No additional storage-off iterations")
    read (on_unit, *, iostat=ios) on_step, on_time, on_iteration, on_norms
    call assert_eq(ios, iostat_end, "No additional storage-on iterations")
    close (off_unit)
    close (on_unit)
  end subroutine compare_history

end program test_residual_storage
