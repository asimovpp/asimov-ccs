!> Verify raw residuals, independent storage and monitoring, and workspace lifetime.
!>
!> On each structured unit-cube mesh, the diagonal matrix is \(A=2I\). The
!> prescribed solution and right-hand side satisfy \(rhs=2\phi+r\), with
!> \(r_{g,j}=j(-1)^g g V_g^{2/3}\) for global cell index \(g\), field component
!> \(j\), and cell volume \(V_g\). Stored residuals must retain these raw signed
!> values, whereas monitoring divides by \(V_g^{2/3}\). For \(N\) global cells,
!> the exact monitored norms are
!> \(RMS_j=j\sqrt{(N+1)(2N+1)/6}\) and \(L_{\infty,j}=jN\).
!>
!> Both matrix evaluation and a precomputed residual exercise four independent
!> storage/monitoring scenarios: neither, monitoring only, storage only, and both.
!> The tests preserve solution and right-hand-side entries, including halos,
!> and distinguish owned-cell norm reductions from halo contributions. Separate
!> checks require snapshots and vector duplicates to own their data and retain
!> their ghost layout. Repeated initialisation and finalisation, followed by a
!> smaller mesh with fewer fields, expose stale workspace and statistics.
program test_residuals

  use testing_lib
  use ccs_base, only: bnd_names_default
  use constants, only: cell
  use core, only: ccs_options, finalise_mesh
  use mesh_utils, only: build_mesh
  use meshing, only: set_mesh_object, nullify_mesh_object, get_local_num_cells, get_total_num_cells, &
                     get_global_num_cells, create_cell_locator, get_global_index, get_volume
  use vec, only: create_vector, duplicate_vector, get_vector_data, restore_vector_data, &
                 get_vector_data_readonly, restore_vector_data_readonly, set_vector_location
  use mat, only: create_matrix, set_nnz, set_matrix_diagonal, finalise_matrix
  use utils, only: initialise, set_size, update
  use residuals, only: compute_residuals, compute_global_residuals, init_residuals, finalise_residuals

  implicit none

  real(ccs_real), parameter :: rtol = 1.0e-12_ccs_real
  real(ccs_real), parameter :: atol = 1.0e-12_ccs_real
  real(ccs_real), parameter :: halo_sentinel = 1.0e12_ccs_real
  real(ccs_real), parameter :: unwritten_residual_sentinel = -1.0e12_ccs_real
  integer(ccs_int), dimension(2), parameter :: mesh_cells_per_side = [4_ccs_int, 3_ccs_int]
  integer(ccs_int), dimension(2), parameter :: case_field_counts = [4_ccs_int, 1_ccs_int]
  character(len=1), dimension(4), parameter :: field_names = ['u', 'v', 'w', 'c']
  type(fluid) :: flow
  type(ccs_residuals) :: residual_statistics
  type(vector_spec) :: vec_properties
  class(ccs_matrix), allocatable :: matrix
  class(ccs_vector), allocatable :: rhs, raw_residual
  integer(ccs_int) :: local_num_cells, total_num_cells, global_num_cells
  integer(ccs_int) :: case_index, field_count

  call init()

  do case_index = 1, size(mesh_cells_per_side)
    field_count = case_field_counts(case_index)
    call initialise_case(mesh_cells_per_side(case_index), field_count)
    call check_duplicate_independence()
    call check_matrix_residuals(store_residuals=.false., monitor_residuals=.false.)
    call check_matrix_residuals(store_residuals=.false., monitor_residuals=.true.)
    call check_matrix_residuals(store_residuals=.true., monitor_residuals=.false.)
    call check_matrix_residuals(store_residuals=.true., monitor_residuals=.true.)
    call check_precomputed_residuals(store_residuals=.false., monitor_residuals=.false.)
    call check_precomputed_residuals(store_residuals=.false., monitor_residuals=.true.)
    call check_precomputed_residuals(store_residuals=.true., monitor_residuals=.false.)
    call check_precomputed_residuals(store_residuals=.true., monitor_residuals=.true.)
    call check_snapshot_independence()
    call finalise_residuals()
    call finalise_residuals()
    call clean_case()
  end do

  call fin()

contains

  !> Construct the diagonal system and independently owned fields on a unit cube.
  !> Reinitialisation must replace the statistics arrays, clear their contents,
  !> and accommodate the reduced field count on the second mesh.
  subroutine initialise_case(cells_per_side, field_count)
    integer(ccs_int), intent(in) :: cells_per_side !< Number of cells along each side of the unit cube.
    integer(ccs_int), intent(in) :: field_count !< Number of velocity/scalar fields in this case.
    type(ccs_options) :: run_options
    type(matrix_spec) :: mat_properties
    real(ccs_real), dimension(:), pointer :: vector_data
    integer(ccs_int) :: field_index

    run_options%mesh%bnd_names = bnd_names_default
    run_options%mesh%cps = cells_per_side
    run_options%mesh%domain_size = 1.0_ccs_real
    mesh = build_mesh(par_env, shared_env, run_options)
    call set_mesh_object(mesh)
    call get_local_num_cells(local_num_cells)
    call get_total_num_cells(total_num_cells)
    call get_global_num_cells(global_num_cells)

    call initialise(vec_properties)
    call set_vector_location(cell, vec_properties)
    call set_size(par_env, mesh, vec_properties)
    call create_vector(vec_properties, rhs)
    call create_vector(vec_properties, raw_residual)
    call get_vector_data(rhs, vector_data)
    vector_data = 2.0_ccs_real
    call restore_vector_data(rhs, vector_data)
    call update(rhs)

    call initialise(mat_properties)
    call set_size(par_env, mesh, mat_properties)
    call set_nnz(1_ccs_int, mat_properties)
    call create_matrix(mat_properties, matrix)
    call set_matrix_diagonal(rhs, matrix)
    call finalise_matrix(matrix)

    allocate (flow%fields(field_count))
    do field_index = 1, field_count
      allocate (field :: flow%fields(field_index)%ptr)
      flow%fields(field_index)%ptr%name = field_names(field_index)
      call create_vector(vec_properties, flow%fields(field_index)%ptr%values)
    end do
    call init_residuals(flow, residual_statistics)
    residual_statistics%L2 = -1.0_ccs_real
    residual_statistics%Linfty = -2.0_ccs_real
    call init_residuals(flow, residual_statistics)
    call assert_eq(size(residual_statistics%L2), field_count, 'Statistics resize with field count')
    call assert_eq(residual_statistics%L2, spread(0.0_ccs_real, 1, field_count), 'Initial RMS statistics reset')
    call assert_eq(residual_statistics%Linfty, spread(0.0_ccs_real, 1, field_count), 'Initial maximum statistics reset')
  end subroutine initialise_case

  !> Check matrix residuals for one independent choice of storage and monitoring.
  !> Fresh output sentinels expose unwritten snapshots; input checks
  !> include every owned and halo entry of both the solution and right-hand side.
  subroutine check_matrix_residuals(store_residuals, monitor_residuals)
    logical, intent(in) :: store_residuals !< Whether each field owns a residual snapshot.
    logical, intent(in) :: monitor_residuals !< Whether the optional statistics argument is supplied.
    integer(ccs_int) :: component

    call configure_storage(store_residuals)
    call init_residuals(flow, residual_statistics)
    do component = 1, field_count
      call set_linear_system(component)
      if (monitor_residuals) then
        call forward_matrix_monitor(flow%fields(component)%ptr, residual_statistics)
      else
        call forward_matrix_monitor(flow%fields(component)%ptr)
      end if
      call check_linear_system_unchanged(component)
      if (store_residuals) call check_stored_residual(component, flow%fields(component)%ptr)
      call assert_bool(allocated(flow%fields(component)%ptr%residuals) .eqv. store_residuals, &
                       'Monitoring preserves requested storage allocation')
    end do
    if (monitor_residuals) call check_global_norms()
  end subroutine check_matrix_residuals

  !> Exercise optional-argument forwarding for the matrix residual interface.
  !> An absent monitor in the caller must remain absent in compute_residuals.
  subroutine forward_matrix_monitor(phi, monitor)
    class(field), intent(inout) :: phi !< Field whose residual may be stored and monitored.
    type(ccs_residuals), optional, intent(inout) :: monitor !< Statistics destination, omitted when monitoring is disabled.

    call compute_residuals(flow, matrix, rhs, phi, monitor)
  end subroutine forward_matrix_monitor

  !> Check precomputed residuals without allowing halo entries into the norms.
  !> Large halo sentinels make an incorrect owned-cell reduction observable. They
  !> remain in the input, while any stored snapshot must refresh its own halos.
  subroutine check_precomputed_residuals(store_residuals, monitor_residuals)
    logical, intent(in) :: store_residuals !< Whether each field owns a residual snapshot.
    logical, intent(in) :: monitor_residuals !< Whether the optional statistics argument is supplied.
    real(ccs_real), dimension(:), pointer :: vector_data
    real(ccs_real), dimension(:), allocatable :: original_residual
    integer(ccs_int) :: index_p, component

    call configure_storage(store_residuals)
    call init_residuals(flow, residual_statistics)
    do component = 1, field_count
      call get_vector_data(raw_residual, vector_data)
      vector_data = halo_sentinel
      do index_p = 1, local_num_cells
        vector_data(index_p) = expected_residual(index_p, component)
      end do
      original_residual = vector_data
      call restore_vector_data(raw_residual, vector_data)
      ! Omit the halo update so the sentinel identifies any halo contribution to a norm.
      if (monitor_residuals) then
        call forward_precomputed_monitor(flow%fields(component)%ptr, residual_statistics)
      else
        call forward_precomputed_monitor(flow%fields(component)%ptr)
      end if
      call get_vector_data_readonly(raw_residual, vector_data)
      call assert_eq(vector_data, original_residual, 'Precomputed residual input including sentinel halos is unchanged')
      call restore_vector_data_readonly(raw_residual, vector_data)
      if (store_residuals) call check_stored_residual(component, flow%fields(component)%ptr)
    end do
    if (monitor_residuals) call check_global_norms()
  end subroutine check_precomputed_residuals

  !> Exercise optional-argument forwarding for the precomputed residual interface.
  !> An absent monitor in the caller must remain absent in compute_residuals.
  subroutine forward_precomputed_monitor(phi, monitor)
    class(field), intent(inout) :: phi !< Field whose residual may be stored and monitored.
    type(ccs_residuals), optional, intent(inout) :: monitor !< Statistics destination, omitted when monitoring is disabled.

    call compute_residuals(flow, raw_residual, phi, monitor)
  end subroutine forward_precomputed_monitor

  !> Require a stored snapshot to survive later use of the private norm workspace.
  !> A separate field with the same name selects the same statistics entry without
  !> owning snapshot storage; two different right-hand sides exercise workspace reuse.
  subroutine check_snapshot_independence()
    class(field), pointer :: phi
    real(ccs_real), dimension(:), pointer :: vector_data

    call configure_storage(.true.)
    phi => flow%fields(1)%ptr
    call set_linear_system(1_ccs_int)
    call compute_residuals(flow, matrix, rhs, phi)
    block
      type(field) :: later_phi
      later_phi%name = phi%name
      call create_vector(vec_properties, later_phi%values)
      call get_vector_data(later_phi%values, vector_data)
      vector_data = 0.0_ccs_real
      call restore_vector_data(later_phi%values, vector_data)
      call update(later_phi%values)
      call compute_residuals(flow, matrix, rhs, later_phi, residual_statistics)
      call check_stored_residual(1_ccs_int, phi)
      ! A changed right-hand side must not overwrite the earlier field snapshot.
      call get_vector_data(rhs, vector_data)
      vector_data = -7.0_ccs_real
      call restore_vector_data(rhs, vector_data)
      call update(rhs)
      call compute_residuals(flow, matrix, rhs, later_phi, residual_statistics)
      call check_stored_residual(1_ccs_int, phi)
    end block
  end subroutine check_snapshot_independence

  !> Verify duplicate ownership, checkout state and the complete ghost layout.
  !> The duplicate is created while the template is checked out. It must then be
  !> independently writable, survive template destruction, and exchange halos
  !> using the original global mapping.
  subroutine check_duplicate_independence()
    class(ccs_vector), allocatable :: template, duplicate
    real(ccs_real), dimension(:), pointer :: vector_data
    integer(ccs_int) :: index_p, global_index_p
    type(cell_locator) :: loc_p

    call create_vector(vec_properties, template)
    call get_vector_data(template, vector_data)
    do index_p = 1, size(vector_data)
      vector_data(index_p) = real(index_p, ccs_real)
    end do
    ! Duplication must not inherit the template's active checkout state.
    call duplicate_vector(template, duplicate)
    call restore_vector_data(template, vector_data)
    call get_vector_data(duplicate, vector_data)
    call assert_eq(size(vector_data), total_num_cells, 'Duplicate preserves ghosted local size')
    vector_data = -3.0_ccs_real
    call restore_vector_data(duplicate, vector_data)
    call get_vector_data_readonly(template, vector_data)
    do index_p = 1, size(vector_data)
      call assert_eq(vector_data(index_p), real(index_p, ccs_real), 'Duplicate does not alias template')
    end do
    call restore_vector_data_readonly(template, vector_data)
    deallocate (template)
    call get_vector_data_readonly(duplicate, vector_data)
    call assert_eq(vector_data, spread(-3.0_ccs_real, 1, total_num_cells), 'Duplicate survives template destruction')
    call restore_vector_data_readonly(duplicate, vector_data)
    call get_vector_data(duplicate, vector_data)
    vector_data = halo_sentinel
    do index_p = 1, local_num_cells
      call create_cell_locator(index_p, loc_p)
      call get_global_index(loc_p, global_index_p)
      vector_data(index_p) = real(global_index_p, ccs_real)
    end do
    call restore_vector_data(duplicate, vector_data)
    call update(duplicate)
    call get_vector_data_readonly(duplicate, vector_data)
    do index_p = 1, total_num_cells
      call create_cell_locator(index_p, loc_p)
      call get_global_index(loc_p, global_index_p)
      call assert_eq(vector_data(index_p), real(global_index_p, ccs_real), 'Duplicate retains ghost mapping')
    end do
    call restore_vector_data_readonly(duplicate, vector_data)
    deallocate (duplicate)
  end subroutine check_duplicate_independence

  !> Recreate the requested storage and initialise every snapshot with a sentinel.
  !> The output sentinel detects entries left unwritten by residual evaluation;
  !> it serves a separate purpose from the large halo values used in norm checks.
  subroutine configure_storage(enabled)
    logical, intent(in) :: enabled !< Whether every field should own residual storage.
    integer(ccs_int) :: field_index
    real(ccs_real), dimension(:), pointer :: vector_data

    do field_index = 1, field_count
      associate (phi => flow%fields(field_index)%ptr)
        if (allocated(phi%residuals)) deallocate (phi%residuals)
        if (enabled) then
          call create_vector(vec_properties, phi%residuals)
          call get_vector_data(phi%residuals, vector_data)
          vector_data = unwritten_residual_sentinel
          call restore_vector_data(phi%residuals, vector_data)
        end if
      end associate
    end do
  end subroutine configure_storage

  !> Set solution and right-hand-side values for the prescribed signed residual.
  !> The component-dependent solution is nonzero, so accidental input changes
  !> cannot be hidden by a trivial solution of the diagonal system.
  subroutine set_linear_system(component)
    integer(ccs_int), intent(in) :: component !< One-based field index setting the residual amplitude.
    real(ccs_real), dimension(:), pointer :: rhs_values, solution_values
    integer(ccs_int) :: index_p

    call get_vector_data(rhs, rhs_values)
    call get_vector_data(flow%fields(component)%ptr%values, solution_values)
    do index_p = 1, total_num_cells
      solution_values(index_p) = expected_solution(index_p, component)
      rhs_values(index_p) = 2.0_ccs_real * solution_values(index_p) &
                            + expected_residual(index_p, component)
    end do
    call restore_vector_data(rhs, rhs_values)
    call restore_vector_data(flow%fields(component)%ptr%values, solution_values)
    call update(rhs)
    call update(flow%fields(component)%ptr%values)
  end subroutine set_linear_system

  !> Check all owned and halo entries of both inputs after residual evaluation.
  subroutine check_linear_system_unchanged(component)
    integer(ccs_int), intent(in) :: component !< One-based field index setting the residual amplitude.
    real(ccs_real), dimension(:), pointer :: vector_data
    integer(ccs_int) :: index_p

    call get_vector_data_readonly(rhs, vector_data)
    do index_p = 1, total_num_cells
      call assert_close(vector_data(index_p), 2.0_ccs_real * expected_solution(index_p, component) &
                        + expected_residual(index_p, component), rtol, atol, 'RHS is unchanged')
    end do
    call restore_vector_data_readonly(rhs, vector_data)
    call get_vector_data_readonly(flow%fields(component)%ptr%values, vector_data)
    do index_p = 1, total_num_cells
      call assert_eq(vector_data(index_p), expected_solution(index_p, component), 'Solution is unchanged')
    end do
    call restore_vector_data_readonly(flow%fields(component)%ptr%values, vector_data)
  end subroutine check_linear_system_unchanged

  !> Verify raw rhs-A phi sign, component values and refreshed halo snapshots.
  subroutine check_stored_residual(component, phi)
    integer(ccs_int), intent(in) :: component !< One-based field index setting the residual amplitude.
    class(field), intent(inout) :: phi !< Field owning the snapshot; inout permits the readonly vector accessor.
    real(ccs_real), dimension(:), pointer :: vector_data
    integer(ccs_int) :: index_p

    call get_vector_data_readonly(phi%residuals, vector_data)
    do index_p = 1, total_num_cells
      call assert_close(vector_data(index_p), expected_residual(index_p, component), rtol, atol, &
                        'Stored raw residual sign, value and halo')
    end do
    call restore_vector_data_readonly(phi%residuals, vector_data)
  end subroutine check_stored_residual

  !> Compare global norms with the analytical sums over owned global cells.
  !> The volume scaling cancels, leaving component times the signed global cell
  !> index. The sum of squares gives the RMS formula stated in the programme header.
  subroutine check_global_norms()
    real(ccs_real) :: expected_rms, global_cell_count
    integer(ccs_int) :: component

    global_cell_count = real(global_num_cells, ccs_real)
    expected_rms = sqrt((global_cell_count + 1.0_ccs_real) * (2.0_ccs_real * global_cell_count + 1.0_ccs_real) / 6.0_ccs_real)
    call compute_global_residuals(par_env, flow, residual_statistics)
    do component = 1, field_count
      call assert_close(residual_statistics%L2(component), real(component, ccs_real) * expected_rms, rtol, atol, &
                        'Global RMS counts owned cells only')
      call assert_close(residual_statistics%Linfty(component), real(component, ccs_real) * global_cell_count, rtol, atol, &
                        'Global maximum counts owned cells only')
    end do
  end subroutine check_global_norms

  !> Evaluate the prescribed signed, volume-scaled residual for an owned or halo cell.
  !> The global index gives the same expectation on each owner and halo copy.
  function expected_residual(index_p, component) result(value)
    integer(ccs_int), intent(in) :: index_p !< Local cell index, including halo entries.
    integer(ccs_int), intent(in) :: component !< One-based field index setting the residual amplitude.
    real(ccs_real) :: value !< Raw signed residual before volume normalisation.
    real(ccs_real) :: volume
    integer(ccs_int) :: global_index_p
    type(cell_locator) :: loc_p

    call create_cell_locator(index_p, loc_p)
    call get_global_index(loc_p, global_index_p)
    call get_volume(loc_p, volume)
    value = real(component * (-1)**global_index_p * global_index_p, ccs_real) * volume**(2.0_ccs_real / 3.0_ccs_real)
  end function expected_residual

  !> Evaluate a nonzero analytical solution for an owned or halo cell.
  !> Distinct component and global-index factors reveal changes to either input.
  function expected_solution(index_p, component) result(value)
    integer(ccs_int), intent(in) :: index_p !< Local cell index, including halo entries.
    integer(ccs_int), intent(in) :: component !< One-based field index setting the residual amplitude.
    real(ccs_real) :: value !< Prescribed solution value for this cell and component.
    integer(ccs_int) :: global_index_p
    type(cell_locator) :: loc_p

    call create_cell_locator(index_p, loc_p)
    call get_global_index(loc_p, global_index_p)
    value = real(global_index_p * (component + 1), ccs_real) / 8.0_ccs_real
  end function expected_solution

  !> Release case vectors and mesh storage in dependency order.
  !> The caller releases private residual workspace first; PETSc vectors are then
  !> destroyed before the mesh shared arrays and the final parallel-environment teardown.
  subroutine clean_case()
    integer(ccs_int) :: field_index

    do field_index = 1, field_count
      deallocate (flow%fields(field_index)%ptr)
    end do
    deallocate (flow%fields, rhs, raw_residual, matrix)
    call nullify_mesh_object()
    call finalise_mesh(par_env, .true.)
  end subroutine clean_case

end program test_residuals
