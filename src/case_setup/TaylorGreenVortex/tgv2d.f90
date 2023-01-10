!v Program file for 2D TaylorGreenVortex case
!
!  @build mpi+petsc

program tgv2d
#include "ccs_macros.inc"

  use petscvec
  use petscsys

  use case_config, only: num_steps, velocity_relax, pressure_relax, res_target, &
                         write_gradients
  use constants, only: cell, face, ccsconfig, ccs_string_len
  use kinds, only: ccs_real, ccs_int
  use types, only: field, upwind_field, central_field, face_field, ccs_mesh, &
                   vector_spec, ccs_vector, field_ptr
  use yaml, only: parse, error_length
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments, sync
  use parallel_types, only: parallel_environment
  use mesh_utils, only: build_square_mesh, write_mesh
  use vec, only: create_vector, set_vector_location
  use petsctypes, only: vector_petsc
  use pv_coupling, only: solve_nonlinear
  use utils, only: set_size, initialise, update, exit_print, calc_kinetic_energy, calc_enstrophy, &
                   add_field_to_outputlist
  use boundary_conditions, only: read_bc_config, allocate_bc_arrays
  use read_config, only: get_bc_variables, get_boundary_count
  use timestepping, only: set_timestep, activate_timestepping, initialise_old_values
  use io_visualisation, only: write_solution
  use fv, only: update_gradient

  implicit none

  class(parallel_environment), allocatable :: par_env
  character(len=:), allocatable :: case_name       ! Case name
  character(len=:), allocatable :: ccs_config_file ! Config file for CCS
  character(len=ccs_string_len), dimension(:), allocatable :: variable_names  ! variable names for BC reading

  type(ccs_mesh) :: mesh
  type(vector_spec) :: vec_properties
  real(ccs_real) :: L

  class(field), allocatable, target :: u, v, w, p, p_prime, mf

  type(field_ptr), allocatable :: output_list(:)

  integer(ccs_int) :: n_boundaries
  integer(ccs_int) :: cps = 50 ! Default value for cells per side

  integer(ccs_int) :: it_start, it_end
  integer(ccs_int) :: irank ! MPI rank ID
  integer(ccs_int) :: isize ! Size of MPI world

  double precision :: start_time
  double precision :: end_time

  logical :: u_sol = .true.  ! Default equations to solve for LDC case
  logical :: v_sol = .true.
  logical :: w_sol = .false.
  logical :: p_sol = .true.

  real(ccs_real) :: dt          ! The timestep
  integer(ccs_int) :: t         ! Timestep counter
  integer(ccs_int) :: nsteps    ! Number of timesteps to perform
  real(ccs_real) :: CFL         ! The CFL target
  integer(ccs_int) :: save_freq ! Frequency of saving solution

  !XXX: Temporary parameters
  integer, parameter :: face_centred = 0         ! Indicates face centred variable
  integer, parameter :: cell_centred_upwind = 1  ! Indicates cell centred variable (upwind scheme)
  integer, parameter :: cell_centred_central = 2 ! Indicates cell centred variable (central scheme)

  ! Launch MPI
  call initialise_parallel_environment(par_env)

  irank = par_env%proc_id
  isize = par_env%num_procs

  call read_command_line_arguments(par_env, cps, case_name=case_name)

  if (irank == par_env%root) print *, "Starting ", case_name, " case!"
  ccs_config_file = case_name // ccsconfig

  call timer(start_time)

  ! Read case name from configuration file
  call read_configuration(ccs_config_file)

  if (irank == par_env%root) then
    call print_configuration()
  end if

  ! Set start and end iteration numbers (eventually will be read from input file)
  it_start = 1
  it_end = num_steps

  ! Create a square mesh
  if (irank == par_env%root) print *, "Building mesh"
  L = 4.0_ccs_real * atan(1.0_ccs_real) ! == PI
  mesh = build_square_mesh(par_env, cps, L)

  ! Initialise fields
  if (irank == par_env%root) print *, "Initialise fields"
  allocate (face_field :: mf)

  ! Create and initialise field vectors
  call initialise(vec_properties)
  call get_boundary_count(ccs_config_file, n_boundaries)
  call get_bc_variables(ccs_config_file, variable_names)

  call set_vector_location(cell, vec_properties)
  call set_size(par_env, mesh, vec_properties)
  call create_field(vec_properties, cell_centred_upwind, "u", u)
  call create_field(vec_properties, cell_centred_upwind, "v", v)
  call create_field(vec_properties, cell_centred_upwind, "w", w)
  call create_field(vec_properties, cell_centred_central, "p", p)
  call create_field(vec_properties, cell_centred_central, "p_prime", p_prime)

  call set_vector_location(face, vec_properties)
  call set_size(par_env, mesh, vec_properties)
  call create_field(vec_properties, 2, "mf", mf)

  ! Add fields to output list
  allocate (output_list(4))
  call add_field_to_outputlist(u, "u", output_list)
  call add_field_to_outputlist(v, "v", output_list)
  call add_field_to_outputlist(w, "w", output_list)
  call add_field_to_outputlist(p, "p", output_list)

  ! Write gradients to solution file
  write_gradients = .true.

  ! Initialise velocity field
  if (irank == par_env%root) print *, "Initialise velocity field"
  call initialise_flow(mesh, u, v, w, p, mf)
  call calc_tgv2d_error(mesh, 0, u, v, w, p)
  call calc_kinetic_energy(par_env, mesh, 0, u, v, w)
  call calc_enstrophy(par_env, mesh, 0, u, v, w)

  ! Solve using SIMPLE algorithm
  if (irank == par_env%root) print *, "Start SIMPLE"

  CFL = 0.1_ccs_real
  dt = CFL * (3.14_ccs_real / 512)
  nsteps = 100
  save_freq = 200

  ! Write out mesh to file
  call write_mesh(par_env, case_name, mesh)

  call activate_timestepping()
  call set_timestep(dt)
  do t = 1, nsteps
    call solve_nonlinear(par_env, mesh, it_start, it_end, res_target, &
                         u_sol, v_sol, w_sol, p_sol, u, v, w, p, p_prime, mf, t)
    call calc_tgv2d_error(mesh, t, u, v, w, p)
    call calc_kinetic_energy(par_env, mesh, t, u, v, w)

    call update_gradient(mesh, u)
    call update_gradient(mesh, v)
    call update_gradient(mesh, w)
    call calc_enstrophy(par_env, mesh, t, u, v, w)

    if ((t == 1) .or. (t == nsteps) .or. (mod(t, save_freq) == 0)) then
      call write_solution(par_env, case_name, mesh, output_list, t, nsteps, dt)
    end if
  end do

  ! Clean-up
  deallocate (u)
  deallocate (v)
  deallocate (w)
  deallocate (p)
  deallocate (p_prime)
  deallocate (output_list)

  call timer(end_time)

  if (irank == par_env%root) then
    print *, "Elapsed time: ", end_time - start_time
  end if

  ! Finalise MPI
  call cleanup_parallel_environment(par_env)

contains

  ! Read YAML configuration file
  subroutine read_configuration(config_filename)

    use read_config, only: get_reference_number, get_steps, &
                           get_convection_scheme, get_relaxation_factor, &
                           get_target_residual

    character(len=*), intent(in) :: config_filename

    class(*), pointer :: config_file_pointer  !< Pointer to CCS config file
    character(len=error_length) :: error

    config_file_pointer => parse(config_filename, error=error)
    if (error /= '') then
      call error_abort(trim(error))
    end if

    call get_steps(config_file_pointer, num_steps)
    if (num_steps == huge(0)) then
      call error_abort("No value assigned to num-steps.")
    end if

    call get_relaxation_factor(config_file_pointer, u_relax=velocity_relax, p_relax=pressure_relax)
    if (velocity_relax == huge(0.0) .and. pressure_relax == huge(0.0)) then
      call error_abort("No values assigned to velocity and pressure underrelaxation.")
    end if

    call get_target_residual(config_file_pointer, res_target)
    if (res_target == huge(0.0)) then
      call error_abort("No value assigned to target residual.")
    end if

  end subroutine

  ! Print test case configuration
  subroutine print_configuration()

    print *, "Solving ", case_name, " case"

    print *, "++++"
    print *, "SIMULATION LENGTH"
    print *, "Running for ", num_steps, "iterations"
    print *, "++++"
    print *, "MESH"
    print *, "Size is ", cps
    print *, "++++"
    print *, "RELAXATION FACTORS"
    write (*, '(1x,a,e10.3)') "velocity: ", velocity_relax
    write (*, '(1x,a,e10.3)') "pressure: ", pressure_relax

  end subroutine

  subroutine initialise_flow(mesh, u, v, w, p, mf)

    use constants, only: insert_mode, ndim
    use types, only: vector_values, cell_locator, face_locator, neighbour_locator
    use meshing, only: set_cell_location, get_global_index, count_neighbours, set_neighbour_location, &
                       get_local_index, set_face_location, get_local_index, get_face_normal, get_centre
    use fv, only: calc_cell_coords
    use utils, only: clear_entries, set_mode, set_row, set_entry, set_values
    use vec, only: get_vector_data, restore_vector_data, create_vector_values

    ! Arguments
    class(ccs_mesh), intent(in) :: mesh
    class(field), intent(inout) :: u, v, w, p, mf

    ! Local variables
    integer(ccs_int) :: n, count
    integer(ccs_int) :: index_p, global_index_p, index_f, index_nb
    real(ccs_real) :: u_val, v_val, w_val, p_val
    type(cell_locator) :: loc_p
    type(face_locator) :: loc_f
    type(neighbour_locator) :: loc_nb
    type(vector_values) :: u_vals, v_vals, w_vals, p_vals
    real(ccs_real), dimension(:), pointer :: mf_data

    real(ccs_real), dimension(ndim) :: x_p, x_f
    real(ccs_real), dimension(ndim) :: face_normal

    integer(ccs_int) :: nnb
    integer(ccs_int) :: j

    ! Set alias
    associate (n_local => mesh%topo%local_num_cells)
      call create_vector_values(n_local, u_vals)
      call create_vector_values(n_local, v_vals)
      call create_vector_values(n_local, w_vals)
      call create_vector_values(n_local, p_vals)
      call set_mode(insert_mode, u_vals)
      call set_mode(insert_mode, v_vals)
      call set_mode(insert_mode, w_vals)
      call set_mode(insert_mode, p_vals)

      ! Set initial values for velocity fields
      do index_p = 1, n_local
        call set_cell_location(mesh, index_p, loc_p)
        call get_global_index(loc_p, global_index_p)

        call get_centre(loc_p, x_p)

        u_val = sin(x_p(1)) * cos(x_p(2))
        v_val = -cos(x_p(1)) * sin(x_p(2))
        w_val = 0.0_ccs_real
        p_val = 0.0_ccs_real !-(sin(2 * x_p(1)) + sin(2 * x_p(2))) * 0.01_ccs_real / 4.0_ccs_real

        call set_row(global_index_p, u_vals)
        call set_entry(u_val, u_vals)
        call set_row(global_index_p, v_vals)
        call set_entry(v_val, v_vals)
        call set_row(global_index_p, w_vals)
        call set_entry(w_val, w_vals)
        call set_row(global_index_p, p_vals)
        call set_entry(p_val, p_vals)
      end do

      call set_values(u_vals, u%values)
      call set_values(v_vals, v%values)
      call set_values(w_vals, w%values)
      call set_values(p_vals, p%values)

      deallocate (u_vals%global_indices)
      deallocate (v_vals%global_indices)
      deallocate (w_vals%global_indices)
      deallocate (p_vals%global_indices)
      deallocate (u_vals%values)
      deallocate (v_vals%values)
      deallocate (w_vals%values)
      deallocate (p_vals%values)
    end associate

    call get_vector_data(mf%values, mf_data)

    count = 0
    n = 0

    ! Loop over local cells and faces
    do index_p = 1, mesh%topo%local_num_cells

      call set_cell_location(mesh, index_p, loc_p)
      call count_neighbours(loc_p, nnb)
      do j = 1, nnb

        call set_neighbour_location(loc_p, j, loc_nb)
        call get_local_index(loc_nb, index_nb)

        ! if neighbour index is greater than previous face index
        if (index_nb > index_p) then ! XXX: abstract this test

          call set_face_location(mesh, index_p, j, loc_f)
          call get_local_index(loc_f, index_f)
          call get_face_normal(loc_f, face_normal)
          call get_centre(loc_f, x_f)

          ! compute initial value based on current face coordinates
          mf_data(index_f) = sin(x_f(1)) * cos(x_f(2)) * face_normal(1) &
                             - cos(x_f(1)) * sin(x_f(2)) * face_normal(2)
        end if

      end do
    end do

    call restore_vector_data(mf%values, mf_data)

    call update(u%values)
    call update(v%values)
    call update(w%values)
    call update(p%values)
    call update(mf%values)

  end subroutine initialise_flow

  subroutine calc_tgv2d_error(mesh, t, u, v, w, p)

    use constants, only: ndim
    use types, only: cell_locator
    use utils, only: str

    use vec, only: get_vector_data, restore_vector_data

    use meshing, only: get_centre, set_cell_location

    use parallel, only: allreduce
    use parallel_types_mpi, only: parallel_environment_mpi

    type(ccs_mesh), intent(in) :: mesh
    integer(ccs_int), intent(in) :: t
    class(field), intent(inout) :: u, v, w, p

    real(ccs_real), dimension(3) :: err_local, err_rms

    real(ccs_real) :: ft
    real(ccs_real) :: u_an, v_an, w_an, p_an
    real(ccs_real), dimension(:), pointer :: u_data, v_data, w_data, p_data

    real(ccs_real) :: mu, rho, nu

    logical, save :: first_time = .true.

    type(cell_locator) :: loc_p
    real(ccs_real), dimension(ndim) :: x_p
    integer(ccs_int) :: index_p

    character(len=ccs_string_len) :: fmt
    real(ccs_real) :: time

    integer :: io_unit

    integer :: ierr

    mu = 0.01_ccs_real ! XXX: currently hardcoded somewhere
    rho = 1.0_ccs_real ! XXX: implicitly 1 throughout
    nu = mu / rho

    err_local(:) = 0.0_ccs_real

    call get_vector_data(u%values, u_data)
    call get_vector_data(v%values, v_data)
    call get_vector_data(w%values, w_data)
    call get_vector_data(p%values, p_data)
    do index_p = 1, mesh%topo%local_num_cells

      call set_cell_location(mesh, index_p, loc_p)
      call get_centre(loc_p, x_p)

      ! Compute analytical solution
      time = t * dt
      ft = exp(-2 * nu * time)
      ! u_an = cos(x_p(1)) * sin(x_p(2)) * ft
      ! v_an = -sin(x_p(1)) * cos(x_p(2)) * ft
      u_an = sin(x_p(1)) * cos(x_p(2)) * ft
      v_an = -cos(x_p(1)) * sin(x_p(2)) * ft
      w_an = 0.0_ccs_real
      p_an = -(rho / 4.0_ccs_real) * (cos(2 * x_p(1)) + cos(2 * x_p(2))) * (ft**2)

      err_local(1) = err_local(1) + (u_an - u_data(index_p))**2
      err_local(2) = err_local(2) + (v_an - v_data(index_p))**2
      err_local(3) = err_local(3) + (w_an - w_data(index_p))**2
      !err_local(4) = err_local(4) + (p_an - p_data(index_p))**2

    end do
    call restore_vector_data(u%values, u_data)
    call restore_vector_data(v%values, v_data)
    call restore_vector_data(w%values, w_data)
    call restore_vector_data(p%values, p_data)

    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_AllReduce(err_local, err_rms, size(err_rms), MPI_DOUBLE, MPI_SUM, par_env%comm, ierr)
    class default
      call error_abort("ERROR: Unknown type")
    end select
    err_rms(:) = sqrt(err_rms(:) / mesh%topo%global_num_cells)

    if (par_env%proc_id == par_env%root) then
      if (first_time) then
        first_time = .false.
        open (newunit=io_unit, file="tgv2d-err.log", status="replace", form="formatted")
      else
        open (newunit=io_unit, file="tgv2d-err.log", status="old", form="formatted", position="append")
      end if
      fmt = '(I0,' // str(size(err_rms)) // '(1x,e12.4))'
      write (io_unit, fmt) t, err_rms
      close (io_unit)
    end if

  end subroutine calc_tgv2d_error

  !v Build a field variable with data and gradient vectors + transient data and boundary arrays.
  subroutine create_field(vec_properties, field_type, field_name, phi)

    use utils, only : debug_print
    
    implicit none
    
    !! Logically vec_properties should be a field_properties variable, but this doesn't yet exist.
    type(vector_spec), intent(in) :: vec_properties !< Vector descriptor for vectors wrapped by field
    integer, intent(in) :: field_type               !< Identifier for what kind of field
    character(len=*), intent(in) :: field_name      !< Field name -- should match against boundary conditions, etc.
    class(field), allocatable, intent(out) :: phi   !< The field being constructed

    call allocate_field(vec_properties, field_type, phi)

    ! XXX: ccs_config_file is host-associated from program scope.
    call read_bc_config(ccs_config_file, field_name, phi)

    !! --- Ensure data is updated/parallel-constructed ---
    ! XXX: Potential abstraction --- see update(vec), etc.
    call update(phi%values)
    if (field_type /= face_centred) then
       ! Current design only computes/stores gradients at cell centres
       call update(phi%x_gradients)
       call update(phi%y_gradients)
       call update(phi%z_gradients)

       call update_gradient(vec_properties%mesh, phi)
    end if
    !! --- End update ---
    
  end subroutine create_field

  !v Allocate a field variable
  subroutine allocate_field(vec_properties, field_type, phi)

    use utils, only : debug_print
    
    implicit none
    
    !! Logically vec_properties should be a field_properties variable, but this doesn't yet exist.
    type(vector_spec), intent(in) :: vec_properties !< Vector descriptor for vectors wrapped by field
    integer, intent(in) :: field_type               !< Identifier for what kind of field
    class(field), allocatable, intent(out) :: phi   !< The field being constructed

    if (field_type == face_centred) then
       call dprint("Create face field")
       allocate (face_field :: phi)
    else if (field_type == cell_centred_upwind) then
       call dprint("Create upwind field")
       allocate (upwind_field :: phi)
    else if (field_type == cell_centred_central) then
       call dprint("Create central field")
       allocate (central_field :: phi)
    end if

    call dprint("Create field values vector")
    call create_vector(vec_properties, phi%values)

    if (field_type /= face_centred) then
       ! Current design only computes/stores gradients at cell centres
       call dprint("Create field gradients vector")
       call create_vector(vec_properties, phi%x_gradients)
       call create_vector(vec_properties, phi%y_gradients)
       call create_vector(vec_properties, phi%z_gradients)

       ! Currently no need for old face values
       call dprint("Create field old values")
       call initialise_old_values(vec_properties, phi)
    end if

    ! XXX: n_boundaries is host-associated from program scope.
    call allocate_bc_arrays(n_boundaries, phi%bcs)

  end subroutine allocate_field
  
end program tgv2d
