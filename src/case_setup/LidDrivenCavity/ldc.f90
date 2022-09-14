!v Program file for LidDrivenCavity case
!
!  @build mpi+petsc

program ldc
#include "ccs_macros.inc"

  use petscvec
  use petscsys

  use case_config, only: num_steps, velocity_relax, pressure_relax, res_target
  use constants, only: cell, face, ccsconfig, ccs_string_len
  use kinds, only: ccs_real, ccs_int
  use types, only: field, upwind_field, central_field, face_field, ccs_mesh, &
                   vector_spec, ccs_vector
  use yaml, only: parse, error_length
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments, sync
  use parallel_types, only: parallel_environment
  use mesh_utils, only: build_square_mesh
  use vec, only: create_vector, set_vector_location
  use petsctypes, only: vector_petsc
  use pv_coupling, only: solve_nonlinear
  use utils, only: set_size, initialise, update, exit_print
  use boundary_conditions, only: read_bc_config, allocate_bc_arrays
  use read_config, only: get_bc_variables, get_boundary_count

  implicit none

  class(parallel_environment), allocatable :: par_env
  character(len=:), allocatable :: case_name       ! Case name
  character(len=:), allocatable :: ccs_config_file ! Config file for CCS
  character(len=ccs_string_len), dimension(:), allocatable :: variable_names  ! variable names for BC reading

  type(ccs_mesh) :: mesh
  type(vector_spec) :: vec_properties

  class(field), allocatable :: u, v, p, p_prime, mf

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

#ifndef EXCLUDE_MISSING_INTERFACE
  integer(ccs_int) :: ierr
  type(tPetscViewer) :: viewer
#endif

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
  mesh = build_square_mesh(par_env, cps, 1.0_ccs_real)

  ! Initialise fields
  if (irank == par_env%root) print *, "Initialise fields"
  allocate (upwind_field :: u)
  allocate (upwind_field :: v)
  allocate (central_field :: p)
  allocate (central_field :: p_prime)
  allocate (face_field :: mf)

  ! Read boundary conditions
  call get_boundary_count(ccs_config_file, n_boundaries)
  call get_bc_variables(ccs_config_file, variable_names)
  call allocate_bc_arrays(n_boundaries, u%bcs)
  call allocate_bc_arrays(n_boundaries, v%bcs)
  call allocate_bc_arrays(n_boundaries, p%bcs)
  call allocate_bc_arrays(n_boundaries, p_prime%bcs)
  call read_bc_config(ccs_config_file, "u", u)
  call read_bc_config(ccs_config_file, "v", v)
  call read_bc_config(ccs_config_file, "p", p)
  call read_bc_config(ccs_config_file, "p", p_prime)

  ! Create and initialise field vectors
  call initialise(vec_properties)

  call set_vector_location(cell, vec_properties)
  call set_size(par_env, mesh, vec_properties)
  call create_vector(vec_properties, u%values)
  call create_vector(vec_properties, v%values)
  call create_vector(vec_properties, p%values)
  call create_vector(vec_properties, p%x_gradients)
  call create_vector(vec_properties, p%y_gradients)
  call create_vector(vec_properties, p%z_gradients)
  call create_vector(vec_properties, p_prime%values)
  call create_vector(vec_properties, p_prime%x_gradients)
  call create_vector(vec_properties, p_prime%y_gradients)
  call create_vector(vec_properties, p_prime%z_gradients)
  call update(u%values)
  call update(v%values)
  call update(p%values)
  call update(p%x_gradients)
  call update(p%y_gradients)
  call update(p%z_gradients)
  call update(p_prime%values)
  call update(p_prime%x_gradients)
  call update(p_prime%y_gradients)
  call update(p_prime%z_gradients)

  call set_vector_location(face, vec_properties)
  call set_size(par_env, mesh, vec_properties)
  call create_vector(vec_properties, mf%values)
  call update(mf%values)

  ! Initialise velocity field
  if (irank == par_env%root) print *, "Initialise velocity field"
  call initialise_velocity(mesh, u, v, mf)
  call update(u%values)
  call update(v%values)
  call update(mf%values)

  ! Solve using SIMPLE algorithm
  if (irank == par_env%root) print *, "Start SIMPLE"
  call solve_nonlinear(par_env, mesh, it_start, it_end, res_target, &
                       u_sol, v_sol, w_sol, p_sol, u, v, p, p_prime, mf)

#ifndef EXCLUDE_MISSING_INTERFACE
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "u", FILE_MODE_WRITE, viewer, ierr)

  associate (vec => u%values)
    select type (vec)
    type is (vector_petsc)
      call VecView(vec%v, viewer, ierr)
    end select
  end associate

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "v", FILE_MODE_WRITE, viewer, ierr)

  associate (vec => v%values)
    select type (vec)
    type is (vector_petsc)
      call VecView(vec%v, viewer, ierr)
    end select
  end associate

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD, "p", FILE_MODE_WRITE, viewer, ierr)

  associate (vec => p%values)
    select type (vec)
    type is (vector_petsc)
      call VecView(vec%v, viewer, ierr)
    end select
  end associate

  call PetscViewerDestroy(viewer, ierr)
#endif

  ! Write out mesh and solution
  call write_solution(par_env, case_name, mesh, cps, u, v, p)

  ! Clean-up
  deallocate (u)
  deallocate (v)
  deallocate (p)
  deallocate (p_prime)

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

  subroutine initialise_velocity(mesh, u, v, mf)

    use constants, only: add_mode
    use types, only: vector_values, cell_locator
    use meshing, only: set_cell_location, get_global_index
    use fv, only: calc_cell_coords
    use utils, only: clear_entries, set_mode, set_row, set_entry, set_values
    use vec, only: get_vector_data, restore_vector_data, create_vector_values

    ! Arguments
    class(ccs_mesh), intent(in) :: mesh
    class(field), intent(inout) :: u, v, mf

    ! Local variables
    integer(ccs_int) :: row, col
    integer(ccs_int) :: index_p, global_index_p
    real(ccs_real) :: u_val, v_val
    type(cell_locator) :: loc_p
    type(vector_values) :: u_vals, v_vals
    real(ccs_real), dimension(:), pointer :: mf_data

    ! Set alias
    associate (n_local => mesh%topo%local_num_cells)
      call create_vector_values(n_local, u_vals)
      call create_vector_values(n_local, v_vals)
      call set_mode(add_mode, u_vals)
      call set_mode(add_mode, v_vals)

      ! Set initial values for velocity fields
      do index_p = 1, n_local
        call set_cell_location(mesh, index_p, loc_p)
        call get_global_index(loc_p, global_index_p)
        call calc_cell_coords(global_index_p, cps, row, col)

        u_val = 0.0_ccs_real
        v_val = 0.0_ccs_real

        call set_row(global_index_p, u_vals)
        call set_entry(u_val, u_vals)
        call set_row(global_index_p, v_vals)
        call set_entry(v_val, v_vals)
      end do

      call set_values(u_vals, u%values)
      call set_values(v_vals, v%values)

      deallocate (u_vals%global_indices)
      deallocate (v_vals%global_indices)
      deallocate (u_vals%values)
      deallocate (v_vals%values)
    end associate

    call get_vector_data(mf%values, mf_data)
    mf_data(:) = 0.0_ccs_real
    call restore_vector_data(mf%values, mf_data)

  end subroutine initialise_velocity

  subroutine write_solution(par_env, case_name, mesh, cps, u, v, p)

    use io, only: initialise_io, cleanup_io, configure_io, &
                  open_file, close_file, write_array
    use kinds, only: ccs_long
    use constants, only: ndim, adiosconfig
    use types, only: io_environment, io_process
    use vec, only : get_vector_data

    ! Arguments
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    character(len=:), allocatable, intent(in) :: case_name
    type(ccs_mesh), intent(in) :: mesh
    integer(ccs_int), intent(in) :: cps
    class(field), intent(in) :: u, v, p

    ! Local variables
    character(len=:), allocatable :: sol_file
    character(len=:), allocatable :: adios2_file
    character(len=:), allocatable :: xdmf_file

    class(io_environment), allocatable :: io_env
    class(io_process), allocatable :: sol_writer

    integer(ccs_long), dimension(1) :: sel_shape
    integer(ccs_long), dimension(1) :: sel_start
    integer(ccs_long), dimension(1) :: sel_count

    integer(ccs_long), dimension(2) :: sel2_shape
    integer(ccs_long), dimension(2) :: sel2_start
    integer(ccs_long), dimension(2) :: sel2_count

    real(ccs_real), dimension(:), pointer :: data

    integer(ccs_int), parameter :: ioxdmf = 999

    sol_file = case_name//'.solution.h5'
    adios2_file = case_name//adiosconfig
    xdmf_file = case_name//'.solution.xmf'

    call initialise_io(par_env, adios2_file, io_env)
    call configure_io(io_env, "sol_writer", sol_writer)

    call open_file(sol_file, "write", sol_writer)

    ! 1D data
    sel_shape(1) = mesh%nglobal
    sel_start(1) = mesh%global_indices(1) - 1
    sel_count(1) = mesh%nlocal

    ! 2D data
    sel2_shape(1) = ndim
    sel2_shape(2) = mesh%nglobal
    sel2_start(1) = 0
    sel2_start(2) = mesh%global_indices(1) - 1
    sel2_count(1) = ndim
    sel2_count(2) = mesh%nlocal

    ! Write mesh cell centre coords
    call write_array(sol_writer, "/xp", sel2_shape, sel2_start, sel2_count, mesh%x_p)

    ! Write u-velocity
    call get_vector_data(u%values, data)
    call write_array(sol_writer, "/u", sel_shape, sel_start, sel_count, data)

    ! Write v-velocity
    call get_vector_data(v%values, data)
    call write_array(sol_writer, "/v", sel_shape, sel_start, sel_count, data)

    ! Write v-velocity
    call get_vector_data(p%values, data)
    call write_array(sol_writer, "/p", sel_shape, sel_start, sel_count, data)

    ! Close the file and ADIOS2 engine
    call close_file(sol_writer)

    ! Finalise the ADIOS2 IO environment
    call cleanup_io(io_env)

    ! Write XML file
    if (par_env%proc_id == par_env%root) then
      ! Open file
      open(ioxdmf, file=xdmf_file, status='unknown')

      ! Write file contents
      write(ioxdmf, '(a)') '<?xml version="1.0"?>'
      write(ioxdmf, '(a)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
      write(ioxdmf, '(a)') '<Xdmf Version="2.0">'
      write(ioxdmf, '(a)') '  <Domain>'
      write(ioxdmf, '(a)') '    <Grid Name="Mesh">'
      write(ioxdmf, '(a,i0,1x,i0,a)') '      <Topology Type="2DSMesh" NumberOfElements="',cps,cps,'"/>'
      write(ioxdmf, '(a)') '      <Geometry Type="XYZ">'
      write(ioxdmf, '(a,i0,1x,i0,1x,i0,a)') '        <DataItem Dimensions="',cps,cps,ndim,'" Format="HDF">'
      write(ioxdmf, '(a,a,a)') '          ',trim(sol_file),':/Step0/xp'
      write(ioxdmf, '(a)') '        </DataItem>'
      write(ioxdmf, '(a)') '      </Geometry>'
      write(ioxdmf, '(a)') '      <Attribute Name="VelocityX" AttributeType="Scalar" Center="Node">'
      write(ioxdmf, '(a,i0,1x,i0,a)') '        <DataItem Dimensions="',cps,cps,'" Format="HDF">'
      write(ioxdmf, '(a,a,a)') '          ',trim(sol_file),':/Step0/u'
      write(ioxdmf, '(a)') '        </DataItem>'
      write(ioxdmf, '(a)') '      </Attribute>'
      write(ioxdmf, '(a)') '      <Attribute Name="VelocityY" AttributeType="Scalar" Center="Node">'
      write(ioxdmf, '(a,i0,1x,i0,a)') '        <DataItem Dimensions="',cps,cps,'" Format="HDF">'
      write(ioxdmf, '(a,a,a)') '          ',trim(sol_file),':/Step0/v'
      write(ioxdmf, '(a)') '        </DataItem>'
      write(ioxdmf, '(a)') '      </Attribute>'
      write(ioxdmf, '(a)') '      <Attribute Name="Pressure" AttributeType="Scalar" Center="Node">'
      write(ioxdmf, '(a,i0,1x,i0,a)') '        <DataItem Dimensions="',cps,cps,'" Format="HDF">'
      write(ioxdmf, '(a,a,a)') '          ',trim(sol_file),':/Step0/p'
      write(ioxdmf, '(a)') '        </DataItem>'
      write(ioxdmf, '(a)') '      </Attribute>'
      write(ioxdmf, '(a)') '    </Grid>'
      write(ioxdmf, '(a)') '  </Domain>'
      write(ioxdmf, '(a)') '</Xdmf>'

      ! Close file
      close(ioxdmf)
    endif ! par_env%root

  end subroutine

end program ldc
