!> @brief Program file for LidDrivenCavity case
!
!> @build mpi+petsc

program ldc

  use petscvec
  use petscsys

  use constants, only : cell, face, ccsconfig
  use bc_constants, only: bc_region_left, bc_region_right, &
                          bc_region_top, bc_region_bottom, &
                          bc_region_live, &
                          bc_type_sym, bc_type_wall
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
  use utils, only: set_size, initialise, update

  implicit none

  class(parallel_environment), allocatable :: par_env
  character(len=:), allocatable :: case_name   !< Case name
  character(len=:), allocatable :: ccs_config_file  !< Config file for CCS

  type(ccs_mesh)    :: square_mesh
  type(vector_spec) :: vec_sizes

  class(field), allocatable :: u, v, p, pp, mf

  integer(ccs_int) :: cps = 50 ! Default value for cells per side

  integer(ccs_int) :: it_start, it_end, ierr
  integer(ccs_int) :: irank !< MPI rank ID
  integer(ccs_int) :: isize !< Size of MPI world

  double precision :: start_time
  double precision :: end_time

  type(tPetscViewer) :: viewer

  ! Reference numbers
  ! real(ccs_real) :: p_ref
  ! real(ccs_real) :: temp_ref
  ! real(ccs_real) :: density
  ! real(ccs_real) :: viscosity
  ! integer(ccs_int) :: pref_at_cell

  ! Number of steps
  integer(ccs_int) :: num_steps

  ! Convection/discretisation scheme
  ! integer(ccs_int) :: u_conv
  ! integer(ccs_int) :: v_conv

  ! Relation factors
  real(ccs_real) :: u_relax
  real(ccs_real) :: p_relax

  ! Boundary conditions - hardcoded for now
  ! integer(ccs_int), dimension(5) :: bnd_region
  ! integer(ccs_int), dimension(5) :: bnd_type
  

  print *, "Starting SIMPLE demo"
  ! Launch MPI
  call initialise_parallel_environment(par_env) 

  irank = par_env%proc_id
  isize = par_env%num_procs

  call read_command_line_arguments(par_env, case_name=case_name)

  ccs_config_file = case_name//ccsconfig

  call timer(start_time)

  ! Read case name from configuration file
  call read_configuration(ccs_config_file)

  ! Set start and end iteration numbers (eventually will be read from input file)
  it_start = 1
  it_end   = num_steps

  ! Create a square mesh
  print *, "Building mesh"
  square_mesh = build_square_mesh(par_env, cps, 1.0_ccs_real)

  ! Initialise fields
  print *, "Initialise fields"
  allocate(upwind_field :: u)
  allocate(upwind_field :: v)
  allocate(central_field :: p)
  allocate(central_field :: pp)
  allocate(face_field :: mf)

  ! Create and initialise field vectors
  call initialise(vec_sizes)

  print *, "Create vectors"
  call set_vector_location(cell, vec_sizes)
  call set_size(par_env, square_mesh, vec_sizes)
  call create_vector(vec_sizes, u%values)
  call create_vector(vec_sizes, v%values)
  call create_vector(vec_sizes, p%values)
  call create_vector(vec_sizes, p%x_gradients)
  call create_vector(vec_sizes, p%y_gradients)
  call create_vector(vec_sizes, p%z_gradients)
  call create_vector(vec_sizes, pp%values)
  call create_vector(vec_sizes, pp%x_gradients)
  call create_vector(vec_sizes, pp%y_gradients)
  call create_vector(vec_sizes, pp%z_gradients)
  call update(u%values)
  call update(v%values)
  call update(p%values)
  call update(p%x_gradients)
  call update(p%y_gradients)
  call update(p%z_gradients)
  call update(pp%values)
  call update(pp%x_gradients)
  call update(pp%y_gradients)
  call update(pp%z_gradients)

  call set_vector_location(face, vec_sizes)
  call set_size(par_env, square_mesh, vec_sizes)
  call create_vector(vec_sizes, mf%values)
  call update(mf%values)
  
  ! Initialise velocity field
  print *, "Initialise velocity field"
  call initialise_velocity(square_mesh, u, v, mf)
  call update(u%values)
  call update(v%values)
  call update(mf%values)

  ! Solve using SIMPLE algorithm
  print *, "Start SIMPLE"
  call solve_nonlinear(par_env, square_mesh, cps, it_start, it_end, u, v, p, pp, mf)

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,"u",FILE_MODE_WRITE,viewer, ierr)

  associate (vec => u%values)
    select type (vec)
    type is (vector_petsc)
      call VecView(vec%v, viewer, ierr)
    end select
  end associate

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,"v",FILE_MODE_WRITE,viewer, ierr)

  associate (vec => v%values)
    select type (vec)
    type is (vector_petsc)
      call VecView(vec%v, viewer, ierr)
    end select
  end associate

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,"p",FILE_MODE_WRITE,viewer, ierr)

  associate (vec => p%values)
    select type (vec)
    type is (vector_petsc)
      call VecView(vec%v, viewer, ierr)
    end select
  end associate

  call PetscViewerDestroy(viewer,ierr)

  ! Clean-up
  deallocate(u)
  deallocate(v)
  deallocate(p)
  deallocate(pp)

  call timer(end_time)

  if(irank == 0) then
     print*, "Elapsed time: ", end_time - start_time
  end if
  
  ! Finalise MPI
  call cleanup_parallel_environment(par_env)

  contains

  ! Read YAML configuration file
  subroutine read_configuration(config_filename)

    use read_config, only: get_reference_number, get_steps, &
                            get_convection_scheme, get_relaxation_factor

    character(len=*), intent(in) :: config_filename
    
    class(*), pointer :: config_file_pointer  !> Pointer to CCS config file
    character(len=error_length) :: error

    config_file_pointer => parse(config_filename, error=error)
    if (error/='') then
      print*,trim(error)
      stop 1
    endif
    
    call get_steps(config_file_pointer, num_steps)
    ! call get_reference_number(config_file_pointer, p_ref=p_ref, temp_ref=temp_ref, &
    !                           dens_ref=density, visc_ref=viscosity, pref_at_cell=pref_at_cell)
    ! call get_convection_scheme(config_file_pointer, u_conv=u_conv, v_conv=v_conv)
    call get_relaxation_factor(config_file_pointer, u_relax=u_relax, p_relax=p_relax)

  end subroutine

  ! Print test case configuration
  ! subroutine print_config()

  !   print*,"Solving ", case_name, " case"

  !   print*,"++++" 
  !   print*,"SIMULATION LENGTH"
  !   print*,"Running for ",num_steps, "time steps"

  !   print*,"++++"   
  !   print*,"REFERENCE NUMBERS"
  !   print*,"Reference pressure ", p_ref, " set a cell ", pref_at_cell
  !   print*,"Reference temperature ", temp_ref
  !   print*,"Density ", density
  !   print*,"Viscocity ", viscosity

  !   print*,"++++" 
  !   print*,"DISCRETISATION SCHEMES"  
  !   if(u_conv == 0) then
  !     print*,"u is upwind" 
  !   else
  !     print*,"u is central" 
  !   end if
  
  !   if(v_conv == 0) then
  !     print*,"v is upwind" 
  !   else
  !     print*,"v is central" 
  !   end if

  !   print*,"++++" 
  !   print*,"RELAXATION FACTORS"
  !   print*,"u is ", u_relax 
  !   print*,"v is ", v_relax 
  !   print*,"p is ", p_relax 

  !   print*,"++++" 
  !   print*,"BOUNDARY CONDITIONS"
  !   print*,"Regions: ", bnd_region
  !   print*,"Types: ", bnd_type

  ! end subroutine

  subroutine initialise_velocity(cell_mesh, u, v, mf)

    use constants, only: add_mode
    use types, only: vector_values, cell_locator
    use meshing, only: set_cell_location, get_global_index
    use fv, only: calc_cell_coords
    use utils, only: pack_entries, set_values
    use vec, only : get_vector_data, restore_vector_data
    
    ! Arguments
    class(ccs_mesh), intent(in) :: cell_mesh
    class(field), intent(inout) :: u, v, mf

    ! Local variables
    integer(ccs_int) :: row, col
    integer(ccs_int) :: local_idx, self_idx
    real(ccs_real) :: u_val, v_val
    type(cell_locator) :: self_loc
    type(vector_values) :: u_vals, v_vals
    real(ccs_real), dimension(:), pointer :: u_data, v_data, mf_data

    ! Set mode
    u_vals%setter_mode = add_mode
    v_vals%setter_mode = add_mode

    ! Set alias
    associate(n_local => cell_mesh%nlocal)
      ! Allocate temporary arrays for storing global cell indices 
      allocate(u_vals%indices(n_local))
      allocate(v_vals%indices(n_local))

      ! Allocate temporary arrays for storing values
      allocate(u_vals%values(n_local))
      allocate(v_vals%values(n_local))

      ! Set initial values for velocity fields
      do local_idx = 1, n_local
        call set_cell_location(cell_mesh, local_idx, self_loc)
        call get_global_index(self_loc, self_idx)
        call calc_cell_coords(self_idx, cps, row, col)

        u_val = real(col, ccs_real)/real(cps, ccs_real)
        v_val = -real(row, ccs_real)/real(cps, ccs_real)

        call pack_entries(local_idx, self_idx, u_val, u_vals)
        call pack_entries(local_idx, self_idx, v_val, v_vals)
      end do
    end associate

    call set_values(u_vals, u%values)
    call set_values(v_vals, v%values)

    deallocate(u_vals%indices)
    deallocate(v_vals%indices)
    deallocate(u_vals%values)
    deallocate(v_vals%values)

    call get_vector_data(u%values, u_data)
    call get_vector_data(v%values, v_data)
    call get_vector_data(mf%values, mf_data)

    u_data(:) = 0.0_ccs_real
    v_data(:) = 0.0_ccs_real
    mf_data(:) = 0.0_ccs_real
    
    call restore_vector_data(u%values, u_data)
    call restore_vector_data(v%values, v_data)
    call restore_vector_data(mf%values, mf_data)
    
  end subroutine initialise_velocity


end program ldc
