program test_distorted_schemes
  use testing_lib

  use bc_constants, only: bc_type_dirichlet
  use constants, only: add_mode, ndim
  use fv, only: calc_diffusion_coeff, calc_mass_flux
  use fv_equations, only: momentum_equation, poisson_equation
  use fv_kernels, only: cd_advection_kernel
  use kinds, only: ccs_int, ccs_real
  use mat, only: create_matrix_values, set_matrix_values_spec_ncols, set_matrix_values_spec_nrows
  use meshing, only: create_cell_locator, create_face_locator, nullify_mesh_object, set_face_interpolation, set_mesh_object
  use types, only: ccs_mesh, central_field, cell_locator, face_locator, matrix_values, matrix_values_spec, vector_values
  use distorted_refs
  use utils, only: set_mode
  use vec, only: create_vector_values

  implicit none

  real(ccs_real), parameter :: rtol = 1.0e-6_ccs_real
  real(ccs_real), parameter :: atol = 1.0e-8_ccs_real

  type(ccs_mesh), target :: mesh_obj
  type(central_field) :: phi
  ! central_field keeps read-only pointers to these arrays during each check.
  real(ccs_real), target :: phi_values(2)
  real(ccs_real), target :: phi_x_gradients(2)
  real(ccs_real), target :: phi_y_gradients(2)
  real(ccs_real), target :: phi_z_gradients(2)
  integer :: icase
  integer :: ifield

  call init()

  call setup_mesh(mesh_obj)
  call setup_field(phi, phi_values, phi_x_gradients, phi_y_gradients, phi_z_gradients)

  do icase = 1, n_cases
    call load_mesh_ref(icase, mesh_obj)
    call check_diffusion(icase)
    do ifield = 1, n_fields
      call load_field_ref(icase, ifield, phi, phi_values, phi_x_gradients, phi_y_gradients, phi_z_gradients)
      call check_mass_flux(icase, ifield)
      call check_momentum(icase, ifield)
      call check_poisson(icase, ifield)
    end do
  end do

  call nullify_mesh_object()
  call fin()

contains

  subroutine check_diffusion(icase)
    integer, intent(in) :: icase

    real(ccs_real) :: coeff
    character(len=32) :: msg

    ! Without cell corrections the coefficient uses the full cell-centre distance.
    call calc_diffusion_coeff(1_ccs_int, 1_ccs_int, .false., mu(1, icase), mu(2, icase), &
                              rho(1, icase), rho(2, icase), schmidt(icase), coeff)
    write (msg, '(a, i0)') "plain diff c", icase
    call assert_close(coeff, diff_plain(icase), rtol, atol, &
                      trim(msg))

    ! With corrections enabled it uses the orthogonal projection onto the face normal.
    call calc_diffusion_coeff(1_ccs_int, 1_ccs_int, .true., mu(1, icase), mu(2, icase), &
                              rho(1, icase), rho(2, icase), schmidt(icase), coeff)
    write (msg, '(a, i0)') "corr diff c", icase
    call assert_close(coeff, diff_corr(icase), rtol, atol, &
                      trim(msg))

  end subroutine check_diffusion

  subroutine check_mass_flux(icase, ifield)
    integer, intent(in) :: icase
    integer, intent(in) :: ifield

    type(face_locator) :: loc_f
    real(ccs_real) :: flux
    character(len=32) :: msg

    call create_face_locator(1_ccs_int, 1_ccs_int, loc_f)
    flux = calc_mass_flux(pressure(:, ifield, icase), grad_p(1, :, ifield, icase), &
                          grad_p(2, :, ifield, icase), grad_p(3, :, ifield, icase), &
                          inv_a(:, icase), loc_f, .true.)

    ! This is the pressure-only Rhie-Chow correction.  It includes the
    ! non-orthogonal pressure-gradient correction and the uneven-cell
    ! interpolation of volume and inverse momentum diagonal.
    write (msg, '(a, i0, a, i0)') "mass c", icase, " f", ifield
    call assert_close(flux, rc_flux(ifield, icase), rtol, atol, &
                      trim(msg))

  end subroutine check_mass_flux

  subroutine check_momentum(icase, ifield)
    integer, intent(in) :: icase
    integer, intent(in) :: ifield

    real(ccs_real) :: zero_flux(1)
    real(ccs_real) :: zero_viscosity(2)
    character(len=32) :: msg

    zero_flux(:) = 0.0_ccs_real
    zero_viscosity(:) = 0.0_ccs_real

    ! These residual checks keep the kernel tests as the coefficient-level
    ! unit tests, while this test verifies equation-object assembly.
    write (msg, '(a, i0, a, i0)') "mom adv c", icase, " f", ifield
    call check_momentum_res(icase, [mflux(icase)], zero_viscosity, &
                            mom_adv_res(ifield, icase), trim(msg))
    write (msg, '(a, i0, a, i0)') "mom diff c", icase, " f", ifield
    call check_momentum_res(icase, zero_flux, mu(:, icase), &
                            mom_diff_res(ifield, icase), trim(msg))
    write (msg, '(a, i0, a, i0)') "mom c", icase, " f", ifield
    call check_momentum_res(icase, [mflux(icase)], mu(:, icase), &
                            mom_res(ifield, icase), trim(msg))

  end subroutine check_momentum

  subroutine check_momentum_res(icase, mass_ref, viscosity_ref, expected, msg)
    integer, intent(in) :: icase
    real(ccs_real), intent(in) :: mass_ref(:)
    real(ccs_real), intent(in) :: viscosity_ref(:)
    real(ccs_real), intent(in) :: expected
    character(len=*), intent(in) :: msg

    type(momentum_equation) :: equation
    type(cd_advection_kernel) :: advection_kernel
    type(matrix_values) :: mat_values
    type(vector_values) :: rhs_values
    type(cell_locator) :: loc_p
    ! momentum_equation stores pointers to these arrays between init and apply.
    real(ccs_real), target :: mass_flux(1)
    real(ccs_real), target :: viscosity(2)
    real(ccs_real), target :: density(2)
    real(ccs_real) :: residual

    mass_flux(:) = mass_ref(:)
    viscosity(:) = viscosity_ref(:)
    density(:) = rho(:, icase)

    call make_work_values(mat_values, rhs_values)
    call equation%init(1_ccs_int, mass_flux, viscosity, density, 1_ccs_int)
    call equation%set_advection(advection_kernel)
    call create_cell_locator(1_ccs_int, loc_p)
    call equation%gather(phi, loc_p)
    call equation%apply(mat_values, rhs_values)

    residual = row_residual(mat_values, rhs_values)
    call assert_close(residual, expected, rtol, atol, msg)

  end subroutine check_momentum_res

  subroutine check_poisson(icase, ifield)
    integer, intent(in) :: icase
    integer, intent(in) :: ifield

    type(poisson_equation) :: equation
    type(matrix_values) :: mat_values
    type(vector_values) :: rhs_values
    type(cell_locator) :: loc_p
    real(ccs_real), target :: inv_a_work(2)
    character(len=32) :: msg

    inv_a_work(:) = inv_a(:, icase)

    call make_work_values(mat_values, rhs_values)
    call equation%init(1_ccs_int)
    call equation%bind_inverse(inv_a_work)
    call create_cell_locator(1_ccs_int, loc_p)
    call equation%gather(phi, loc_p)
    call equation%apply(mat_values, rhs_values)

    ! The pressure-correction equation is diffusion-like; check the assembled
    ! residual rather than repeating the diffusion-kernel coefficient checks.
    write (msg, '(a, i0, a, i0)') "pois c", icase, " f", ifield
    call assert_close(row_residual(mat_values, rhs_values), pc_res(ifield, icase), rtol, atol, &
                      trim(msg))

  end subroutine check_poisson

  subroutine setup_mesh(mesh_obj)
    type(ccs_mesh), target, intent(inout) :: mesh_obj

    call set_mesh_object(mesh_obj)

    mesh_obj%topo%global_num_cells = 2_ccs_int
    mesh_obj%topo%local_num_cells = 2_ccs_int
    mesh_obj%topo%halo_num_cells = 0_ccs_int
    mesh_obj%topo%total_num_cells = 2_ccs_int
    mesh_obj%topo%global_num_faces = 1_ccs_int
    mesh_obj%topo%num_faces = 1_ccs_int
    mesh_obj%topo%max_faces = 1_ccs_int
    mesh_obj%topo%shared_array_local_offset = 0_ccs_int
    mesh_obj%topo%shared_array_total_offset = 0_ccs_int
    mesh_obj%is_generated = .true.
    mesh_obj%geo%h = 1.0_ccs_real
    mesh_obj%geo%scalefactor = 1.0_ccs_real

    allocate (mesh_obj%topo%global_indices(2))
    allocate (mesh_obj%topo%natural_indices(2))
    allocate (mesh_obj%topo%global_face_indices(1, 2))
    allocate (mesh_obj%topo%face_indices(1, 2))
    allocate (mesh_obj%topo%nb_indices(1, 2))
    allocate (mesh_obj%topo%num_nb(2))
    allocate (mesh_obj%topo%face_cell1(1))
    allocate (mesh_obj%topo%face_cell2(1))
    allocate (mesh_obj%topo%bnd_rid(1))

    mesh_obj%topo%global_indices(:) = [1_ccs_int, 2_ccs_int]
    mesh_obj%topo%natural_indices(:) = [1_ccs_int, 2_ccs_int]
    mesh_obj%topo%global_face_indices(:, :) = 1_ccs_int
    mesh_obj%topo%face_indices(:, :) = 1_ccs_int
    mesh_obj%topo%nb_indices(1, :) = [2_ccs_int, 1_ccs_int]
    mesh_obj%topo%num_nb(:) = 1_ccs_int
    mesh_obj%topo%face_cell1(1) = 1_ccs_int
    mesh_obj%topo%face_cell2(1) = 2_ccs_int
    mesh_obj%topo%bnd_rid(1) = 0_ccs_int

    allocate (mesh_obj%geo%x_p(ndim, 2))
    allocate (mesh_obj%geo%x_f(ndim, 1, 2))
    allocate (mesh_obj%geo%face_normals(ndim, 1, 2))
    allocate (mesh_obj%geo%face_areas(1, 2))
    allocate (mesh_obj%geo%volumes(2))
    allocate (mesh_obj%geo%face_interpol(1))

  end subroutine setup_mesh

  subroutine load_mesh_ref(icase, mesh_obj)
    integer, intent(in) :: icase
    type(ccs_mesh), target, intent(inout) :: mesh_obj

    type(face_locator) :: loc_f

    mesh_obj%geo%h = dx_mag(icase)
    mesh_obj%geo%x_p(:, 1) = xp(:, icase)
    mesh_obj%geo%x_p(:, 2) = xn(:, icase)
    mesh_obj%geo%x_f(:, 1, 1) = xf(:, icase)
    mesh_obj%geo%x_f(:, 1, 2) = xf(:, icase)
    mesh_obj%geo%face_normals(:, 1, 1) = normal(:, icase)
    mesh_obj%geo%face_normals(:, 1, 2) = -normal(:, icase)
    mesh_obj%geo%face_areas(:, :) = area(icase)
    mesh_obj%geo%volumes(:) = volume(:, icase)
    mesh_obj%geo%face_interpol(:) = 0.0_ccs_real

    call create_face_locator(1_ccs_int, 1_ccs_int, loc_f)
    call set_face_interpolation(lf(icase), loc_f)

  end subroutine load_mesh_ref

  subroutine setup_field(phi, values, x_gradients, y_gradients, z_gradients)
    type(central_field), intent(inout) :: phi
    real(ccs_real), target, intent(inout) :: values(:)
    real(ccs_real), target, intent(inout) :: x_gradients(:)
    real(ccs_real), target, intent(inout) :: y_gradients(:)
    real(ccs_real), target, intent(inout) :: z_gradients(:)

    phi%values_ro => values
    phi%x_gradients_ro => x_gradients
    phi%y_gradients_ro => y_gradients
    phi%z_gradients_ro => z_gradients
    phi%enable_cell_corrections = .true.

    allocate (phi%bcs%bc_types(1))
    phi%bcs%bc_types(1) = bc_type_dirichlet

  end subroutine setup_field

  subroutine load_field_ref(icase, ifield, phi, values, x_gradients, y_gradients, z_gradients)
    integer, intent(in) :: icase
    integer, intent(in) :: ifield
    type(central_field), intent(inout) :: phi
    real(ccs_real), intent(inout) :: values(:)
    real(ccs_real), intent(inout) :: x_gradients(:)
    real(ccs_real), intent(inout) :: y_gradients(:)
    real(ccs_real), intent(inout) :: z_gradients(:)

    values(:) = phi_val(:, ifield, icase)
    x_gradients(:) = grad_phi(1, :, ifield, icase)
    y_gradients(:) = grad_phi(2, :, ifield, icase)
    z_gradients(:) = grad_phi(3, :, ifield, icase)
    phi%Schmidt = schmidt(icase)

  end subroutine load_field_ref

  subroutine make_work_values(mat_values, rhs_values)
    type(matrix_values), intent(out) :: mat_values
    type(vector_values), intent(out) :: rhs_values

    type(matrix_values_spec) :: mat_spec

    call set_matrix_values_spec_nrows(1_ccs_int, mat_spec)
    call set_matrix_values_spec_ncols(2_ccs_int, mat_spec)
    call create_matrix_values(mat_spec, mat_values)
    call create_vector_values(1_ccs_int, rhs_values)
    call set_mode(add_mode, mat_values)
    call set_mode(add_mode, rhs_values)

  end subroutine make_work_values

  real(ccs_real) function row_residual(mat_values, rhs_values) result(residual)
    type(matrix_values), intent(in) :: mat_values
    type(vector_values), intent(in) :: rhs_values

    residual = rhs_values%values(1)
    residual = residual - matrix_entry(mat_values, 1_ccs_int) * phi_values(1)
    residual = residual - matrix_entry(mat_values, 2_ccs_int) * phi_values(2)

  end function row_residual

  real(ccs_real) function matrix_entry(mat_values, global_col) result(entry)
    type(matrix_values), intent(in) :: mat_values
    integer(ccs_int), intent(in) :: global_col

    integer(ccs_int), dimension(1) :: col_index

    col_index = findloc(mat_values%global_col_indices, global_col - 1_ccs_int, kind=ccs_int)
    if (col_index(1) == 0) then
      error stop "Requested matrix column was not assembled"
    end if

    entry = mat_values%values(col_index(1))

  end function matrix_entry

end program test_distorted_schemes
