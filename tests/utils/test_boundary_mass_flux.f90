!> Verify reusable boundary mass-flux integration in serial and across MPI ranks.
program test_boundary_mass_flux

  use core, only: ccs_options
  use testing_lib
  use flow_stats, only: boundary_flux_context, boundary_mass_flux_balance, &
                        create_boundary_flux_context, &
                        integrate_boundary_mass_fluxes, report_boundary_mass_flux_balance
  use kinds, only: ccs_int, ccs_real
  use logging, only: log_unit_out
  use meshing, only: nullify_mesh_object, set_mesh_object
  use timestepping, only: activate_timestepping, finalise_timestep, &
                          reset_timestepping, set_timestep
  use types, only: ccs_mesh, face_field, patch_type_generic, &
                   patch_type_inflow, patch_type_outflow

  implicit none

  real(ccs_real), parameter :: atol = 0.0_ccs_real
  real(ccs_real), parameter :: rtol = 0.0_ccs_real
  character(len=32) :: test_mode

  call init()

  call get_command_argument(1, test_mode)
  if (len_trim(test_mode) > 0) then
    if (par_env%num_procs /= 1) then
      call stop_test("Boundary mass-flux validation tests require one MPI rank")
    end if
    call test_validation(trim(test_mode))
    call stop_test("Expected boundary mass-flux validation failure")
  end if

  select case (par_env%num_procs)
  case (1)
    call test_serial_integration()
    call test_context_reconstruction()
    call test_zero_patches()
    call test_reporting()
    call test_empty_reporting()
  case (4)
    call test_four_rank_integration()
    call test_four_rank_no_boundary_faces()
    call test_four_rank_zero_patches()
  case default
    call stop_test("Boundary mass-flux test requires one or four MPI ranks")
  end select

  call fin()

contains

  subroutine test_validation(mode)
    character(len=*), intent(in) :: mode

    type(boundary_flux_context) :: context
    type(boundary_mass_flux_balance) :: balance
    type(ccs_mesh) :: test_mesh
    type(face_field) :: mass_flux_f
    real(ccs_real), dimension(:), allocatable, target :: values

    call build_local_mesh(test_mesh, 1_ccs_int, [1_ccs_int], [2.0_ccs_real])

    select case (mode)
    case ("uninitialized")
      values = [3.0_ccs_real]
      mass_flux_f%values_ro => values
      balance = integrate_boundary_mass_fluxes(par_env, mass_flux_f, context)
    case ("disassociated-values")
      context = create_test_boundary_flux_context(test_mesh)
      balance = integrate_boundary_mass_fluxes(par_env, mass_flux_f, context)
    case ("face-extent")
      allocate (values(0))
      mass_flux_f%values_ro => values
      context = create_test_boundary_flux_context(test_mesh)
      balance = integrate_boundary_mass_fluxes(par_env, mass_flux_f, context)
    case default
      call stop_test("Unknown boundary mass-flux validation mode")
    end select
  end subroutine test_validation

  subroutine test_serial_integration()
    type(boundary_flux_context) :: context
    type(boundary_mass_flux_balance) :: balance
    type(ccs_mesh) :: test_mesh
    type(face_field) :: mass_flux_f
    real(ccs_real), dimension(:), allocatable, target :: values

    call build_local_mesh(test_mesh, 3_ccs_int, &
                          [1_ccs_int, 2_ccs_int, 1_ccs_int, 3_ccs_int], &
                          [2.0_ccs_real, 4.0_ccs_real, 0.5_ccs_real, 3.0_ccs_real])
    values = [1.5_ccs_real, 4.0_ccs_real, &
              (0.0_ccs_real - 1.0_ccs_real), 3.0_ccs_real]
    mass_flux_f%values_ro => values

    context = create_test_boundary_flux_context(test_mesh)
    balance = integrate_boundary_mass_fluxes(par_env, mass_flux_f, context)

    call assert_eq(size(balance%patches), 3, "Expected all boundary patches")
    call assert_eq(balance%patches(1)%name, "patch-1", "Expected patch metadata")
    call assert_close(balance%patches%integrated_mass_flux, &
                      [8.0_ccs_real, (0.0_ccs_real - 4.0_ccs_real), 4.5_ccs_real], &
                      rtol, atol, "Expected signed area-weighted patch totals")
    call assert_close(balance%global_imbalance, 8.5_ccs_real, rtol, atol, &
                      "Expected nonzero global imbalance")

    values = [(0.0_ccs_real - 1.0_ccs_real), 6.0_ccs_real, &
              0.5_ccs_real, (0.0_ccs_real - 2.0_ccs_real)]
    balance = integrate_boundary_mass_fluxes(par_env, mass_flux_f, context)
    call assert_close(balance%patches%integrated_mass_flux, &
                      [(0.0_ccs_real - 1.0_ccs_real), 2.0_ccs_real, &
                       (0.0_ccs_real - 3.0_ccs_real)], &
                      rtol, atol, "Expected context reuse with updated face values")
    call assert_close(balance%global_imbalance, 0.0_ccs_real - 2.0_ccs_real, &
                      rtol, atol, &
                      "Expected reused-context imbalance")

    nullify (mass_flux_f%values_ro)
  end subroutine test_serial_integration

  subroutine test_context_reconstruction()
    type(boundary_flux_context) :: context
    type(boundary_mass_flux_balance) :: balance
    type(ccs_mesh) :: test_mesh
    type(face_field) :: mass_flux_f
    real(ccs_real), dimension(:), allocatable, target :: values

    call build_local_mesh(test_mesh, 1_ccs_int, [1_ccs_int], [2.0_ccs_real])
    values = [3.0_ccs_real]
    mass_flux_f%values_ro => values
    context = create_test_boundary_flux_context(test_mesh)

    call build_local_mesh(test_mesh, 2_ccs_int, [2_ccs_int, 1_ccs_int], &
                          [3.0_ccs_real, 4.0_ccs_real])
    test_mesh%boundary_patches(1)%name = "reconstructed-inlet"
    test_mesh%boundary_patches(1)%patch_type = patch_type_inflow
    test_mesh%boundary_patches(2)%name = "reconstructed-outlet"
    test_mesh%boundary_patches(2)%patch_type = patch_type_outflow
    nullify (mass_flux_f%values_ro)
    values = [5.0_ccs_real, (0.0_ccs_real - 2.0_ccs_real)]
    mass_flux_f%values_ro => values
    context = create_test_boundary_flux_context(test_mesh)
    balance = integrate_boundary_mass_fluxes(par_env, mass_flux_f, context)

    call assert_eq(size(balance%patches), 2, "Expected reconstructed patch count")
    call assert_eq(balance%patches(1)%name, "reconstructed-inlet", &
                   "Expected reconstructed patch metadata")
    call assert_eq(balance%patches(2)%patch_type, patch_type_outflow, &
                   "Expected reconstructed patch classification")
    call assert_close(balance%patches%integrated_mass_flux, &
                      [20.0_ccs_real, (0.0_ccs_real - 6.0_ccs_real)], &
                      rtol, atol, "Expected reconstructed topology and geometry")
    call assert_close(balance%global_imbalance, 14.0_ccs_real, rtol, atol, &
                      "Expected reconstructed-context imbalance")

    nullify (mass_flux_f%values_ro)
  end subroutine test_context_reconstruction

  subroutine test_zero_patches()
    type(boundary_flux_context) :: context
    type(boundary_mass_flux_balance) :: balance
    type(ccs_mesh) :: test_mesh
    type(face_field) :: mass_flux_f
    integer(ccs_int), dimension(:), allocatable :: no_patch_ids
    real(ccs_real), dimension(:), allocatable :: no_areas
    real(ccs_real), dimension(:), allocatable, target :: values

    allocate (no_patch_ids(0), no_areas(0))
    call build_local_mesh(test_mesh, 0_ccs_int, no_patch_ids, no_areas)
    allocate (values(0))
    mass_flux_f%values_ro => values

    context = create_test_boundary_flux_context(test_mesh)
    balance = integrate_boundary_mass_fluxes(par_env, mass_flux_f, context)

    call assert_eq(size(balance%patches), 0, "Expected an empty patch result")
    call assert_close(balance%global_imbalance, 0.0_ccs_real, rtol, atol, &
                      "Expected zero imbalance for an empty boundary")

    nullify (mass_flux_f%values_ro)
  end subroutine test_zero_patches

  subroutine test_four_rank_integration()
    type(boundary_flux_context) :: context
    type(boundary_mass_flux_balance) :: balance
    type(ccs_mesh) :: test_mesh
    type(face_field) :: mass_flux_f
    integer(ccs_int), dimension(:), allocatable :: no_patch_ids
    real(ccs_real), dimension(:), allocatable :: no_areas
    real(ccs_real), dimension(:), allocatable, target :: values

    allocate (no_patch_ids(0), no_areas(0))
    select case (par_env%proc_id)
    case (0)
      call build_local_mesh(test_mesh, 4_ccs_int, no_patch_ids, no_areas)
      allocate (values(0))
    case (1)
      call build_local_mesh(test_mesh, 4_ccs_int, [1_ccs_int], [1.5_ccs_real])
      values = [2.0_ccs_real]
    case (2)
      call build_local_mesh(test_mesh, 4_ccs_int, [1_ccs_int, 2_ccs_int], &
                            [2.0_ccs_real, 0.5_ccs_real])
      values = [4.0_ccs_real, (0.0_ccs_real - 1.0_ccs_real)]
    case (3)
      call build_local_mesh(test_mesh, 4_ccs_int, [3_ccs_int], [2.0_ccs_real])
      values = [(0.0_ccs_real - 3.0_ccs_real)]
    end select
    mass_flux_f%values_ro => values

    context = create_test_boundary_flux_context(test_mesh)
    balance = integrate_boundary_mass_fluxes(par_env, mass_flux_f, context)

    call assert_eq(size(balance%patches), 4, "Expected every patch on every rank")
    call assert_close(balance%patches%integrated_mass_flux, &
                      [1.0_ccs_real, 2.0_ccs_real, &
                       (0.0_ccs_real - 6.0_ccs_real), 0.0_ccs_real], &
                      rtol, atol, "Expected dense patch reduction")
    call assert_close(balance%global_imbalance, 0.0_ccs_real - 3.0_ccs_real, &
                      rtol, atol, &
                      "Expected identical imbalance on every rank")

    nullify (mass_flux_f%values_ro)
  end subroutine test_four_rank_integration

  subroutine test_four_rank_no_boundary_faces()
    type(boundary_flux_context) :: context
    type(boundary_mass_flux_balance) :: balance
    type(ccs_mesh) :: test_mesh
    type(face_field) :: mass_flux_f
    integer(ccs_int), dimension(:), allocatable :: no_patch_ids
    real(ccs_real), dimension(:), allocatable :: no_areas
    real(ccs_real), dimension(:), allocatable, target :: values

    allocate (no_patch_ids(0), no_areas(0), values(0))
    call build_local_mesh(test_mesh, 2_ccs_int, no_patch_ids, no_areas)
    mass_flux_f%values_ro => values

    context = create_test_boundary_flux_context(test_mesh)
    balance = integrate_boundary_mass_fluxes(par_env, mass_flux_f, context)

    call assert_eq(size(balance%patches), 2, "Expected configured empty patches")
    call assert_close(balance%patches%integrated_mass_flux, &
                      [0.0_ccs_real, 0.0_ccs_real], rtol, atol, &
                      "Expected zero totals without owned boundary faces")
    call assert_close(balance%global_imbalance, 0.0_ccs_real, rtol, atol, &
                      "Expected zero all-rank empty-boundary imbalance")

    nullify (mass_flux_f%values_ro)
  end subroutine test_four_rank_no_boundary_faces

  subroutine test_four_rank_zero_patches()
    type(boundary_flux_context) :: context
    type(boundary_mass_flux_balance) :: balance
    type(ccs_mesh) :: test_mesh
    type(face_field) :: mass_flux_f
    integer(ccs_int), dimension(:), allocatable :: no_patch_ids
    real(ccs_real), dimension(:), allocatable :: no_areas
    real(ccs_real), dimension(:), allocatable, target :: values

    allocate (no_patch_ids(0), no_areas(0), values(0))
    call build_local_mesh(test_mesh, 0_ccs_int, no_patch_ids, no_areas)
    mass_flux_f%values_ro => values

    context = create_test_boundary_flux_context(test_mesh)
    balance = integrate_boundary_mass_fluxes(par_env, mass_flux_f, context)

    call assert_eq(size(balance%patches), 0, "Expected no configured patches")
    call assert_close(balance%global_imbalance, 0.0_ccs_real, rtol, atol, &
                      "Expected zero all-rank zero-patch imbalance")

    nullify (mass_flux_f%values_ro)
  end subroutine test_four_rank_zero_patches

  subroutine test_reporting()
    type(boundary_mass_flux_balance) :: balance
    type(ccs_options) :: run_options
    integer :: i
    integer :: report_unit
    integer :: saved_log_unit

    allocate (balance%patches(2))
    balance%patches(1)%name = "boundary-name-that-must-not-be-truncated"
    balance%patches(1)%patch_type = patch_type_outflow
    balance%patches(1)%integrated_mass_flux = 1.25_ccs_real
    balance%patches(2)%name = "inlet"
    balance%patches(2)%patch_type = patch_type_inflow
    balance%patches(2)%integrated_mass_flux = 0.0_ccs_real - 2.5_ccs_real
    balance%global_imbalance = 0.0_ccs_real - 1.25_ccs_real
    run_options%diagnostics%boundary_mass_fluxes = .true.

    call reset_timestepping()
    call activate_timestepping()
    call set_timestep(0.05_ccs_real)
    do i = 1, 7
      call finalise_timestep()
    end do

    saved_log_unit = log_unit_out
    open (newunit=report_unit, status="scratch", action="readwrite")
    log_unit_out = report_unit
    call report_boundary_mass_flux_balance(par_env, run_options, balance)
    rewind (report_unit)

    call assert_report_line(report_unit, &
                            "Boundary mass-flux balance: step 7, time 3.500000E-01")
    call assert_report_line(report_unit, &
                            "Sign convention: positive outward, negative inward")
    call assert_report_line(report_unit, repeat("-", 94))
    call assert_report_line(report_unit, &
                            "Boundary ID | Boundary name" // repeat(" ", 27) //  &
                            " | Patch type | Mass flux [solver units]")
    call assert_report_line(report_unit, repeat("-", 94))
    call assert_report_line(report_unit, &
                            "1" // repeat(" ", 10) //  &
                            " | boundary-name-that-must-not-be-truncated" //  &
                            " | outflow" // repeat(" ", 3) //  &
                            " | " // repeat(" ", 11) // "+1.250000E+00")
    call assert_report_line(report_unit, &
                            "2" // repeat(" ", 10) //  &
                            " | inlet" // repeat(" ", 35) //  &
                            " | inflow" // repeat(" ", 4) //  &
                            " | " // repeat(" ", 11) // "-2.500000E+00")
    call assert_report_line(report_unit, repeat("-", 94))
    call assert_report_line(report_unit, &
                            repeat(" ", 11) //  &
                            " | GLOBAL IMBALANCE" // repeat(" ", 24) //  &
                            " | " // repeat(" ", 10) //  &
                            " | " // repeat(" ", 11) // "-1.250000E+00")
    call assert_report_line(report_unit, repeat("-", 94))
    call assert_end_of_report(report_unit)

    close (report_unit)
    log_unit_out = saved_log_unit
    call reset_timestepping()
  end subroutine test_reporting

  subroutine test_empty_reporting()
    type(boundary_mass_flux_balance) :: balance
    type(ccs_options) :: run_options
    integer :: report_unit
    integer :: saved_log_unit

    allocate (balance%patches(0))
    balance%global_imbalance = 0.0_ccs_real
    run_options%diagnostics%boundary_mass_fluxes = .true.

    saved_log_unit = log_unit_out
    open (newunit=report_unit, status="scratch", action="readwrite")
    log_unit_out = report_unit
    call report_boundary_mass_flux_balance(par_env, run_options, balance)
    rewind (report_unit)

    call assert_report_line(report_unit, &
                            "Boundary mass-flux balance: steady solve")
    call assert_report_line(report_unit, &
                            "Sign convention: positive outward, negative inward")
    call assert_report_line(report_unit, repeat("-", 70))
    call assert_report_line(report_unit, &
                            "Boundary ID | Boundary name    | Patch type" //  &
                            " | Mass flux [solver units]")
    call assert_report_line(report_unit, repeat("-", 70))
    call assert_report_line(report_unit, "(no boundary patches)")
    call assert_report_line(report_unit, repeat("-", 70))
    call assert_report_line(report_unit, &
                            repeat(" ", 11) //  &
                            " | GLOBAL IMBALANCE | " // repeat(" ", 10) //  &
                            " | " // repeat(" ", 11) // "+0.000000E+00")
    call assert_report_line(report_unit, repeat("-", 70))
    call assert_end_of_report(report_unit)

    close (report_unit)
    log_unit_out = saved_log_unit
  end subroutine test_empty_reporting

  subroutine assert_report_line(unit, expected)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: expected

    character(len=512) :: actual
    integer :: stat

    read (unit, '(a)', iostat=stat) actual
    if (stat /= 0) call stop_test("Expected another boundary mass-flux report line")
    if (trim(actual) /= expected) then
      print *, "Expected report line: [", expected, "]"
      print *, "Actual report line:   [", trim(actual), "]"
      call stop_test("Boundary mass-flux report line did not match")
    end if
  end subroutine assert_report_line

  subroutine assert_end_of_report(unit)
    integer, intent(in) :: unit

    character(len=512) :: unexpected
    integer :: stat

    read (unit, '(a)', iostat=stat) unexpected
    if (stat == 0) call stop_test("Boundary mass-flux report had unexpected extra lines")
  end subroutine assert_end_of_report

  function create_test_boundary_flux_context(test_mesh) result(context)
    type(ccs_mesh), target, intent(inout) :: test_mesh
    type(boundary_flux_context) :: context

    call set_mesh_object(test_mesh)
    context = create_boundary_flux_context(test_mesh)
    call nullify_mesh_object()
  end function create_test_boundary_flux_context

  subroutine build_local_mesh(test_mesh, patch_count, patch_ids, areas)
    type(ccs_mesh), intent(out) :: test_mesh
    integer(ccs_int), intent(in) :: patch_count
    integer(ccs_int), dimension(:), intent(in) :: patch_ids
    real(ccs_real), dimension(:), intent(in) :: areas

    integer(ccs_int) :: face_count
    integer(ccs_int) :: ghost_cell
    integer(ccs_int) :: i
    integer(ccs_int) :: max_faces

    face_count = size(patch_ids, kind=ccs_int)
    call assert_eq(size(areas), face_count, "Synthetic face metadata must agree")

    allocate (test_mesh%boundary_patches(patch_count))
    do i = 1, patch_count
      write (test_mesh%boundary_patches(i)%name, '(a,i0)') "patch-", i
      test_mesh%boundary_patches(i)%patch_type = patch_type_generic
    end do

    test_mesh%topo%local_num_cells = merge(1_ccs_int, 0_ccs_int, face_count > 0)
    test_mesh%topo%total_num_cells = test_mesh%topo%local_num_cells + &
                                     merge(1_ccs_int, 0_ccs_int, patch_count > 0)
    test_mesh%topo%num_faces = face_count
    max_faces = max(1_ccs_int, face_count)
    test_mesh%topo%max_faces = max_faces
    test_mesh%topo%shared_array_local_offset = 1_ccs_int

    allocate (test_mesh%topo%num_nb(test_mesh%topo%total_num_cells), source=0_ccs_int)
    allocate (test_mesh%topo%nb_indices(max_faces, test_mesh%topo%total_num_cells), &
              source=0_ccs_int)
    allocate (test_mesh%topo%face_indices(max_faces, test_mesh%topo%total_num_cells), &
              source=0_ccs_int)
    allocate (test_mesh%geo%face_areas( &
              max_faces, test_mesh%topo%total_num_cells + &
              test_mesh%topo%shared_array_local_offset), source=101.0_ccs_real)

    if (face_count > 0) then
      test_mesh%topo%num_nb(1) = face_count
      test_mesh%topo%nb_indices(:, 1) = -patch_ids
      test_mesh%topo%face_indices(:, 1) = [(face_count - i + 1, i=1, face_count)]
      test_mesh%geo%face_areas(:, 1 + test_mesh%topo%shared_array_local_offset) = areas
    end if

    if (patch_count > 0) then
      ghost_cell = test_mesh%topo%total_num_cells
      test_mesh%topo%num_nb(ghost_cell) = 1_ccs_int
      test_mesh%topo%nb_indices(1, ghost_cell) = -1_ccs_int
      test_mesh%topo%face_indices(1, ghost_cell) = max(1_ccs_int, face_count)
      test_mesh%geo%face_areas(1, ghost_cell + &
                               test_mesh%topo%shared_array_local_offset) = 17.0_ccs_real
    end if
  end subroutine build_local_mesh

end program test_boundary_mass_flux
